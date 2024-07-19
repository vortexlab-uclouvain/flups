/**
 * Implementation of functions that are common to different GPU backends
 * Mostly concerns the real-to-real transforms.
 *
 */
#include "backend.hpp"

#if M_HAVE_GPU

flups_fft_plan flups_plan_r2r(flups_fft_arch arch, double *in, double *out, flups_kind kind) {
    flups_fft_plan plan;
    // parameters for the r2r transform
    plan.n       = arch.dim;
    plan.dist    = arch.dist;
    plan.stride  = arch.stride;
    plan.howmany = arch.howmany;
    plan.kind    = kind;
    plan.is_r2r  = true;
    // check whether we need to plan a r2c or c2r (by default it's a r2c)
    bool is_r2c = true;

    switch (kind) {
        // DCT Type I
        case FFTW_REDFT00:
            plan.dim     = 2 * plan.n - 2;
            plan.idist   = 2 * plan.n;  // should be 2n-2 (but cufft doesn't like it)
            plan.odist   = plan.n;
            plan.memsize = (2 * plan.n) * plan.howmany;
            break;
        // DST Type I
        case FFTW_RODFT00:
            plan.dim     = 2 * plan.n + 2;
            plan.idist   = 2 * plan.n + 4;  // should be 2n+2 (but cufft doesn't like it)
            plan.odist   = plan.n + 2;
            plan.memsize = (2 * plan.n + 4) * plan.howmany;
            break;
        // DCT/DST Type II
        case FFTW_REDFT10:
        case FFTW_RODFT10:
            plan.dim     = 2 * plan.n;
            plan.idist   = 2 * plan.n + 2;  // should be 2n (but cufft doesn't like it)
            plan.odist   = plan.n + 1;
            plan.memsize = (2 * plan.n + 2) * plan.howmany;
            break;
        // DCT/DST Type III
        case FFTW_REDFT01:
        case FFTW_RODFT01:
            plan.dim     = 2 * plan.n;
            plan.idist   = plan.n + 1;
            plan.odist   = 2 * plan.n + 2;  // here too !
            plan.memsize = (2 * plan.n + 2) * plan.howmany;
            // this is a c2r transform
            is_r2c = false;
            break;
        // DCT/DST Type IV
        case FFTW_REDFT11:
        case FFTW_RODFT11:
            plan.dim     = 8 * plan.n;
            plan.idist   = 8 * plan.n + 2;  //  should be 8n (but cufft doesn't like it)
            plan.odist   = 4 * plan.n + 1;
            plan.memsize = (8 * plan.n + 2) * plan.howmany;
            break;
        default:
            FLUPS_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }

    // Note: some plans (typically involved in the computation of the Green function
    //       for mixed BC's) require a non-unit stride of 2. Here, we handle that
    //       directly during the pre- and post-processing (thus the plan has unit stride).
    m_cufft_check_errors(cufftPlanMany(&plan.fft, 1, &plan.dim,
                                       &plan.dim, 1, plan.idist,
                                       &plan.dim, 1, plan.odist,
                                       is_r2c ? CUFFT_D2Z : CUFFT_Z2D, plan.howmany));
    return plan;
}

void flups_execute_r2r(flups_fft_plan plan, double *in, double *out) {
    int n       = plan.n;
    int dim     = plan.dim;
    int dist    = plan.dist;
    int idist   = plan.idist;
    int odist   = plan.odist;
    int stride  = plan.stride;
    int howmany = plan.howmany;

    // Working memory
    double *mem = plan.mem;

    switch (plan.kind) {
        // DCT Type I
        case FFTW_REDFT00: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + i] = in[j * dist + stride * i];

                if (i == 0 || i == n - 1) {
                    mem[j * idist + i] *= .5;
                }
            };
            for_loop(pre_op, n * howmany, "DCT Type I pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = 2.0 * mem[2 * (j * odist + i)];
            };
            for_loop(post_op, n * howmany, "DCT Type I post-processing");

            break;
        }
        // DCT Type II
        case FFTW_REDFT10: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + i] = in[j * dist + stride * i];
            };
            for_loop(pre_op, n * howmany, "DCT Type II pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                double t                   = -i * M_PI / (2.0 * n);
                out[j * dist + stride * i] = 2.0 * (mem[2 * (j * odist + i)] * cos(t) - mem[2 * (j * odist + i) + 1] * sin(t));
            };
            for_loop(post_op, n * howmany, "DCT Type II post-processing");

            break;
        }
        // DCT Type III
        case FFTW_REDFT01: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                double t                     = -i * M_PI / (2.0 * n);
                mem[2 * (j * idist + i)]     = in[j * dist + stride * i] * cos(t);
                mem[2 * (j * idist + i) + 1] = in[j * dist + stride * i] * sin(t);
            };
            for_loop(pre_op, n * howmany, "DCT Type III pre-processing");

            flups_execute_c2r(plan, (flups_complex *)mem, mem);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = mem[1 + j * odist + i];
            };
            for_loop(post_op, n * howmany, "DCT Type III post-processing");

            break;
        }
        // DCT Type IV
        case FFTW_REDFT11: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + (2 * i + 1)]           = in[j * dist + stride * i];
                mem[j * idist + (4 * n - (2 * i + 1))] = -in[j * dist + stride * i];
            };
            for_loop(pre_op, n * howmany, "DCT Type IV pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = mem[2 * (j * odist + (2 * i + 1))];
            };
            for_loop(post_op, n * howmany, "DCT Type IV post-processing");

            break;
        }
        // DST Type I
        case FFTW_RODFT00: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[1 + j * idist + i] = in[j * dist + stride * i];
            };
            for_loop(pre_op, n * howmany, "DST Type I pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = -2.0 * mem[3 + 2 * (j * odist + i)];
            };
            for_loop(post_op, n * howmany, "DST Type I post-processing");

            break;
        }
        // DST Type II
        case FFTW_RODFT10: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                const int sign     = 2 * (i % 2 == 0) - 1;
                mem[j * idist + i] = sign * in[j * dist + stride * i];
            };
            for_loop(pre_op, n * howmany, "DST Type II pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                double t                             = -i * M_PI / (2.0 * n);
                out[j * dist + (n - 1) - stride * i] = 2.0 * (mem[2 * (j * odist + i)] * cos(t) - mem[2 * (j * odist + i) + 1] * sin(t));
            };
            for_loop(post_op, n * howmany, "DST Type II post-processing");

            break;
        }
        // DST Type III
        case FFTW_RODFT01: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                double t                     = -i * M_PI / (2.0 * n);
                mem[2 * (j * idist + i)]     = in[j * dist + (n - 1) - stride * i] * cos(t);
                mem[2 * (j * idist + i) + 1] = in[j * dist + (n - 1) - stride * i] * sin(t);
            };
            for_loop(pre_op, n * howmany, "DST Type III pre-processing");

            flups_execute_c2r(plan, (flups_complex *)mem, mem);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                const int sign             = 2 * (i % 2 == 0) - 1;
                out[j * dist + stride * i] = sign * mem[1 + j * odist + i];
            };
            for_loop(post_op, n * howmany, "DST Type III post-processing");

            break;
        }
        // DST Type IV
        case FFTW_RODFT11: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + (2 * i + 1)]           = in[j * dist + (n - 1) - stride * i];
                mem[j * idist + (4 * n - (2 * i + 1))] = -in[j * dist + (n - 1) - stride * i];
            };
            for_loop(pre_op, n * howmany, "DST Type IV pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                const int sign             = 2 * (i % 2 == 0) - 1;
                out[j * dist + stride * i] = sign * mem[2 * (j * odist + (2 * i + 1))];
            };
            for_loop(post_op, n * howmany, "DST Type IV post-processing");

            break;
        }
        default:
            FLUPS_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }

    // Reset memory
    backend_memset(mem, 0, plan.memsize * sizeof(double));
}

#endif  // M_HAVE_GPU
