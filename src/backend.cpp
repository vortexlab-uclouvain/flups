/**
 * Implementation of functions that are common to different GPU backends
 * Mostly concerns the real-to-real transforms.
 *
 */
#include "backend.hpp"

#if M_HAVE_GPU

flups_fft_plan flups_plan_r2r(flups_fft_arch arch, double *in, double *out, flups_r2r_kind_t kind) {
    flups_fft_plan plan;
    // parameters for the r2r transform
    plan.n       = arch.dim;
    plan.dist    = arch.dist;
    plan.stride  = arch.stride;
    plan.howmany = arch.howmany;
    plan.kind    = kind;
    plan.is_r2r  = true;
    // check which underlying transform we need to plan (by default it's a r2c)
    enum { XFORM_R2C,
           XFORM_C2R,
           XFORM_C2C } xform = XFORM_R2C;

    switch (kind) {
        // DCT Type I
        case FLUPS_DCT_I:
            plan.dim     = 2 * plan.n - 2;
            plan.idist   = 2 * plan.n;  // 2*(dim/2+1): in-place r2c padding, dim=2n-2
            plan.odist   = plan.n;
            plan.memsize = (2 * plan.n) * plan.howmany;
            break;
        // DST Type I
        case FLUPS_DST_I:
            plan.dim     = 2 * plan.n + 2;
            plan.idist   = 2 * plan.n + 4;  // 2*(dim/2+1): in-place r2c padding, dim=2n+2
            plan.odist   = plan.n + 2;
            plan.memsize = (2 * plan.n + 4) * plan.howmany;
            break;
        // DCT/DST Type II: Makhoul's fast algorithm, a single n-point real FFT
        // (even/odd-reordered input + complex post-twiddle) instead of a 2n-point
        // zero-padded one -- halves the scratch memory vs. the naive embedding.
        case FLUPS_DCT_II:
        case FLUPS_DST_II:
            plan.dim     = plan.n;
            plan.idist   = plan.n + 2;      // padded for in-place r2c (real units)
            plan.odist   = plan.n / 2 + 1;  // complex units
            plan.memsize = (plan.n + 2) * plan.howmany;
            break;
        // DCT/DST Type III: adjoint of Type II, a single n-point real IFFT
        // (complex pre-twiddle + de-reordering) instead of a 2n-point one.
        case FLUPS_DCT_III:
        case FLUPS_DST_III:
            plan.dim     = plan.n;
            plan.idist   = plan.n / 2 + 1;  // complex units
            plan.odist   = plan.n + 2;      // padded for in-place c2r (real units)
            plan.memsize = (plan.n + 2) * plan.howmany;
            // this is a c2r transform
            xform = XFORM_C2R;
            break;
        // DCT/DST Type IV, n even: even/odd-split x into a,b (n/2 samples each),
        // pack them as one complex signal (a + i*b), pre-twiddle, and run a single
        // n-point complex FFT; A[k] and B[k] (the transforms of a and b alone) are
        // then recovered from C[k] and C[n-1-k] via a conjugate-symmetry relation
        // that holds because a and b are separately real. This halves the scratch
        // memory again vs. the plain 2n-point twiddle-FFT approach below
        // (8n -> 4n -> 2n reals). n odd falls back to the 2n-point version, since
        // the even/odd split does not divide evenly.
        case FLUPS_DCT_IV:
        case FLUPS_DST_IV:
            if (plan.n % 2 == 0) {
                plan.dim     = plan.n;
                plan.idist   = plan.n;                       // complex units, exact
                plan.odist   = plan.n;                       // complex units
                plan.memsize = (2 * plan.n) * plan.howmany;  // n complex numbers = 2n doubles
            } else {
                plan.dim     = 2 * plan.n;
                plan.idist   = 2 * plan.n;                   // complex units, exact
                plan.odist   = 2 * plan.n;                   // complex units
                plan.memsize = (4 * plan.n) * plan.howmany;  // 2n complex numbers = 4n doubles
            }
            xform = XFORM_C2C;
            break;
        default:
            BACKEND_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }

        // Note: some plans (typically involved in the computation of the Green function
        //       for mixed BC's) require a non-unit stride of 2. Here, we handle that
        //       directly during the pre- and post-processing (thus the plan has unit stride).
#if defined(BACKEND_CUDA)
    cufftType cuda_type = (xform == XFORM_R2C) ? CUFFT_D2Z : (xform == XFORM_C2R) ? CUFFT_Z2D
                                                                                  : CUFFT_Z2Z;
    m_cufft_check_errors(cufftPlanMany(&plan.fft, 1, &plan.dim,
                                       &plan.dim, 1, plan.idist,
                                       &plan.dim, 1, plan.odist,
                                       cuda_type, plan.howmany));
#elif defined(BACKEND_HIP)
    hipfftType hip_type = (xform == XFORM_R2C) ? HIPFFT_D2Z : (xform == XFORM_C2R) ? HIPFFT_Z2D
                                                                                   : HIPFFT_Z2Z;
    m_hipfft_check_errors(hipfftPlanMany(&plan.fft, 1, &plan.dim,
                                         &plan.dim, 1, plan.idist,
                                         &plan.dim, 1, plan.odist,
                                         hip_type, plan.howmany));
#else
#error "Unsupported backend, something went wrong"
#endif
    return plan;
}

void flups_execute_r2r(flups_fft_plan plan, double *in, double *out, stream_t stream) {
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
        case FLUPS_DCT_I: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + i] = in[j * dist + stride * i];

                if (i == 0 || i == n - 1) {
                    mem[j * idist + i] *= .5;
                }
            };
            parallel_for(pre_op, n * howmany, "DCT Type I pre-processing", stream);

            flups_execute_r2c(plan, mem, (flups_complex *)mem, stream);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = 2.0 * mem[2 * (j * odist + i)];
            };
            parallel_for(post_op, n * howmany, "DCT Type I post-processing", stream);

            break;
        }
        // DCT Type II -- Makhoul's algorithm: even/odd-reorder the n samples into
        // v, take a single n-point real FFT of v, then extract X[k] from V[min(k,n-k)]
        // with a complex post-twiddle (using conjugate symmetry for k > n/2).
        case FLUPS_DCT_II: {
            // Pre-processing: v[p] = x[2p] for p < ceil(n/2), v[p] = x[2(n-1-p)+1] otherwise
            auto pre_op = m_backend_lambda(const int idx) {
                const int p    = idx % n;
                const int j    = idx / n;
                const int half = (n + 1) / 2;
                const int src  = (p < half) ? (2 * p) : (2 * (n - 1 - p) + 1);

                mem[j * idist + p] = in[j * dist + stride * src];
            };
            parallel_for(pre_op, n * howmany, "DCT Type II pre-processing", stream);

            flups_execute_r2c(plan, mem, (flups_complex *)mem, stream);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int    k    = idx % n;
                const int    j    = idx / n;
                const int    half = n / 2;
                const int    m    = (k <= half) ? k : (n - k);
                const double sign = (k <= half) ? 1.0 : -1.0;
                const double t    = M_PI * k / (2.0 * n);
                const double re   = mem[2 * (j * odist + m)];
                const double im   = sign * mem[2 * (j * odist + m) + 1];

                out[j * dist + stride * k] = 2.0 * (re * cos(t) + im * sin(t));
            };
            parallel_for(post_op, n * howmany, "DCT Type II post-processing", stream);

            break;
        }
        // DCT Type III -- adjoint of Type II: rebuild the complex spectrum V from
        // pairs (X[k], X[n-k]) via the (self-inverse) twiddle rotation, take a single
        // n-point real IFFT to get v, then de-reorder v back into x.
        case FLUPS_DCT_III: {
            const int half = n / 2;
            const int hp1  = half + 1;

            // Pre-processing: rebuild V[m], m = 0..n/2
            auto pre_op = m_backend_lambda(const int idx) {
                const int m = idx % hp1;
                const int j = idx / hp1;

                double re, im;
                if (m == 0) {
                    re = 0.5 * in[j * dist];
                    im = 0.0;
                } else {
                    const int    kc  = n - m;
                    const double t   = M_PI * m / (2.0 * n);
                    const double Xm  = in[j * dist + stride * m];
                    const double Xkc = in[j * dist + stride * kc];
                    re               = 0.5 * (cos(t) * Xm + sin(t) * Xkc);
                    im               = 0.5 * (sin(t) * Xm - cos(t) * Xkc);
                }
                mem[2 * (j * idist + m)]     = re;
                mem[2 * (j * idist + m) + 1] = im;
            };
            parallel_for(pre_op, hp1 * howmany, "DCT Type III pre-processing", stream);

            flups_execute_c2r(plan, (flups_complex *)mem, mem, stream);

            // Post-processing: de-reorder v (a factor 2 restores FFTW's REDFT01
            // convention: REDFT10 followed by REDFT01 reproduces the input scaled by 2n)
            auto post_op = m_backend_lambda(const int idx) {
                const int p    = idx % n;
                const int j    = idx / n;
                const int half = (n + 1) / 2;
                const int dst  = (p < half) ? (2 * p) : (2 * (n - 1 - p) + 1);

                out[j * dist + stride * dst] = 2.0 * mem[j * odist + p];
            };
            parallel_for(post_op, n * howmany, "DCT Type III post-processing", stream);

            break;
        }
        // DCT Type IV
        case FLUPS_DCT_IV: {
            if (n % 2 == 0) {
                // n even: even/odd-split x into a[m]=x[2m], b[m]=x[2m+1] (m=0..n/2-1),
                // pack c[m]=a[m]+i*b[m], pre-twiddle, single n-point complex FFT.
                // Post-processing recovers A[k]=DFT-like(a), B[k]=DFT-like(b) from
                // C[k] and C[n-1-k] (valid because a, b are separately real), then
                // combines them with the usual DCT-IV twiddle.
                const int M      = n / 2;
                auto      pre_op = m_backend_lambda(const int idx) {
                    const int    m = idx % M;
                    const int    j = idx / M;
                    const double a = in[j * dist + stride * (2 * m)];
                    const double b = in[j * dist + stride * (2 * m + 1)];
                    const double t = M_PI * m / n;

                    mem[2 * (j * idist + m)]     = a * cos(t) - b * sin(t);
                    mem[2 * (j * idist + m) + 1] = a * sin(t) + b * cos(t);
                };
                parallel_for(pre_op, M * howmany, "DCT Type IV pre-processing", stream);

                flups_execute_c2c(plan, (flups_complex *)mem, (flups_complex *)mem, 1, stream);

                auto post_op = m_backend_lambda(const int idx) {
                    const int    k     = idx % n;
                    const int    j     = idx / n;
                    const double Cr    = mem[2 * (j * odist + k)];
                    const double Ci    = mem[2 * (j * odist + k) + 1];
                    const double C2r   = mem[2 * (j * odist + (n - 1 - k))];
                    const double C2i   = mem[2 * (j * odist + (n - 1 - k)) + 1];
                    const double Ca    = 0.5 * (Cr + C2r);
                    const double Sa    = 0.5 * (Ci - C2i);
                    const double Cb    = 0.5 * (Ci + C2i);
                    const double Sb    = 0.5 * (C2r - Cr);
                    const double theta = M_PI * (2 * k + 1) / (4.0 * n);

                    out[j * dist + stride * k] = 2.0 * (Ca * cos(theta) - Sa * sin(theta) +
                                                        Cb * cos(3.0 * theta) - Sb * sin(3.0 * theta));
                };
                parallel_for(post_op, n * howmany, "DCT Type IV post-processing", stream);
            } else {
                // n odd fallback: complex pre-twiddle w[i]=x[i]*exp(i*pi*i/2n), a single
                // 2n-point complex FFT (inverse-direction kernel, sign=+1) of the
                // n-nonzero, n-zero-padded w, then a complex post-twiddle.
                auto pre_op = m_backend_lambda(const int idx) {
                    const int    i  = idx % n;
                    const int    j  = idx / n;
                    const double t  = M_PI * i / (2.0 * n);
                    const double xi = in[j * dist + stride * i];

                    mem[2 * (j * idist + i)]     = xi * cos(t);
                    mem[2 * (j * idist + i) + 1] = xi * sin(t);
                };
                parallel_for(pre_op, n * howmany, "DCT Type IV pre-processing", stream);

                flups_execute_c2c(plan, (flups_complex *)mem, (flups_complex *)mem, 1, stream);

                auto post_op = m_backend_lambda(const int idx) {
                    const int    k   = idx % n;
                    const int    j   = idx / n;
                    const double phi = M_PI * (2 * k + 1) / (4.0 * n);
                    const double re  = mem[2 * (j * odist + k)];
                    const double im  = mem[2 * (j * odist + k) + 1];

                    out[j * dist + stride * k] = 2.0 * (re * cos(phi) - im * sin(phi));
                };
                parallel_for(post_op, n * howmany, "DCT Type IV post-processing", stream);
            }

            break;
        }
        // DST Type I
        case FLUPS_DST_I: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[1 + j * idist + i] = in[j * dist + stride * i];
            };
            parallel_for(pre_op, n * howmany, "DST Type I pre-processing", stream);

            flups_execute_r2c(plan, mem, (flups_complex *)mem, stream);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = -2.0 * mem[3 + 2 * (j * odist + i)];
            };
            parallel_for(post_op, n * howmany, "DST Type I post-processing", stream);

            break;
        }
        // DST Type II -- DST-II(x)[k] = DCT-II(x')[n-1-k] with x'_i = (-1)^i x_i:
        // same reordering as DCT-II with a sign flip on the source sample, and the
        // twiddle-extracted value written to the mirrored output position.
        case FLUPS_DST_II: {
            // Pre-processing
            auto pre_op = m_backend_lambda(const int idx) {
                const int    p    = idx % n;
                const int    j    = idx / n;
                const int    half = (n + 1) / 2;
                const int    src  = (p < half) ? (2 * p) : (2 * (n - 1 - p) + 1);
                const double sign = (src % 2 == 0) ? 1.0 : -1.0;

                mem[j * idist + p] = sign * in[j * dist + stride * src];
            };
            parallel_for(pre_op, n * howmany, "DST Type II pre-processing", stream);

            flups_execute_r2c(plan, mem, (flups_complex *)mem, stream);

            // Post-processing
            auto post_op = m_backend_lambda(const int idx) {
                const int    k    = idx % n;
                const int    j    = idx / n;
                const int    half = n / 2;
                const int    m    = (k <= half) ? k : (n - k);
                const double sign = (k <= half) ? 1.0 : -1.0;
                const double t    = M_PI * k / (2.0 * n);
                const double re   = mem[2 * (j * odist + m)];
                const double im   = sign * mem[2 * (j * odist + m) + 1];

                out[j * dist + (n - 1) - stride * k] = 2.0 * (re * cos(t) + im * sin(t));
            };
            parallel_for(post_op, n * howmany, "DST Type II post-processing", stream);

            break;
        }
        // DST Type III -- DST-III(X)[i] = (-1)^i * DCT-III(X')[i] with X'_k = X[n-1-k]:
        // same spectrum rebuild as DCT-III fed with the mirrored input, sign-flipped
        // on the final (de-reordered) spatial output.
        case FLUPS_DST_III: {
            const int half = n / 2;
            const int hp1  = half + 1;

            // Pre-processing: rebuild V[m] from the mirrored input X'_k = X[n-1-k]
            auto pre_op = m_backend_lambda(const int idx) {
                const int m = idx % hp1;
                const int j = idx / hp1;

                double re, im;
                if (m == 0) {
                    re = 0.5 * in[j * dist + (n - 1)];
                    im = 0.0;
                } else {
                    const int    kc  = n - m;
                    const double t   = M_PI * m / (2.0 * n);
                    const double Xm  = in[j * dist + (n - 1) - stride * m];
                    const double Xkc = in[j * dist + (n - 1) - stride * kc];
                    re               = 0.5 * (cos(t) * Xm + sin(t) * Xkc);
                    im               = 0.5 * (sin(t) * Xm - cos(t) * Xkc);
                }
                mem[2 * (j * idist + m)]     = re;
                mem[2 * (j * idist + m) + 1] = im;
            };
            parallel_for(pre_op, hp1 * howmany, "DST Type III pre-processing", stream);

            flups_execute_c2r(plan, (flups_complex *)mem, mem, stream);

            // Post-processing: de-reorder v and apply the (-1)^i sign (see DCT-III for the factor 2)
            auto post_op = m_backend_lambda(const int idx) {
                const int    p    = idx % n;
                const int    j    = idx / n;
                const int    half = (n + 1) / 2;
                const int    dst  = (p < half) ? (2 * p) : (2 * (n - 1 - p) + 1);
                const double sign = (dst % 2 == 0) ? 1.0 : -1.0;

                out[j * dist + stride * dst] = 2.0 * sign * mem[j * odist + p];
            };
            parallel_for(post_op, n * howmany, "DST Type III post-processing", stream);

            break;
        }
        // DST Type IV -- DST-IV(x)[k] = (-1)^k * DCT-IV(x_reversed)[k], x_reversed[i] = x[n-1-i]:
        // same pre/post-twiddle machinery as DCT-IV fed with the reversed input, with a
        // (-1)^k sign flip on the output (same index, no position reversal here).
        case FLUPS_DST_IV: {
            if (n % 2 == 0) {
                // a[m] = x_reversed[2m] = x[n-1-2m], b[m] = x_reversed[2m+1] = x[n-2-2m]
                const int M      = n / 2;
                auto      pre_op = m_backend_lambda(const int idx) {
                    const int    m = idx % M;
                    const int    j = idx / M;
                    const double a = in[j * dist + (n - 1) - stride * (2 * m)];
                    const double b = in[j * dist + (n - 1) - stride * (2 * m + 1)];
                    const double t = M_PI * m / n;

                    mem[2 * (j * idist + m)]     = a * cos(t) - b * sin(t);
                    mem[2 * (j * idist + m) + 1] = a * sin(t) + b * cos(t);
                };
                parallel_for(pre_op, M * howmany, "DST Type IV pre-processing", stream);

                flups_execute_c2c(plan, (flups_complex *)mem, (flups_complex *)mem, 1, stream);

                auto post_op = m_backend_lambda(const int idx) {
                    const int    k     = idx % n;
                    const int    j     = idx / n;
                    const double Cr    = mem[2 * (j * odist + k)];
                    const double Ci    = mem[2 * (j * odist + k) + 1];
                    const double C2r   = mem[2 * (j * odist + (n - 1 - k))];
                    const double C2i   = mem[2 * (j * odist + (n - 1 - k)) + 1];
                    const double Ca    = 0.5 * (Cr + C2r);
                    const double Sa    = 0.5 * (Ci - C2i);
                    const double Cb    = 0.5 * (Ci + C2i);
                    const double Sb    = 0.5 * (C2r - Cr);
                    const double theta = M_PI * (2 * k + 1) / (4.0 * n);
                    const double sign  = (k % 2 == 0) ? 1.0 : -1.0;

                    out[j * dist + stride * k] = sign * 2.0 * (Ca * cos(theta) - Sa * sin(theta) + Cb * cos(3.0 * theta) - Sb * sin(3.0 * theta));
                };
                parallel_for(post_op, n * howmany, "DST Type IV post-processing", stream);
            } else {
                auto pre_op = m_backend_lambda(const int idx) {
                    const int    i  = idx % n;
                    const int    j  = idx / n;
                    const double t  = M_PI * i / (2.0 * n);
                    const double xi = in[j * dist + (n - 1) - stride * i];

                    mem[2 * (j * idist + i)]     = xi * cos(t);
                    mem[2 * (j * idist + i) + 1] = xi * sin(t);
                };
                parallel_for(pre_op, n * howmany, "DST Type IV pre-processing", stream);

                flups_execute_c2c(plan, (flups_complex *)mem, (flups_complex *)mem, 1, stream);

                auto post_op = m_backend_lambda(const int idx) {
                    const int    k    = idx % n;
                    const int    j    = idx / n;
                    const double phi  = M_PI * (2 * k + 1) / (4.0 * n);
                    const double re   = mem[2 * (j * odist + k)];
                    const double im   = mem[2 * (j * odist + k) + 1];
                    const double sign = (k % 2 == 0) ? 1.0 : -1.0;

                    out[j * dist + stride * k] = sign * 2.0 * (re * cos(phi) - im * sin(phi));
                };
                parallel_for(post_op, n * howmany, "DST Type IV post-processing", stream);
            }

            break;
        }
        default:
            BACKEND_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }

    // Reset memory
    backend_memset(mem, 0, plan.memsize * sizeof(double), stream);
}

#endif  // M_HAVE_GPU
