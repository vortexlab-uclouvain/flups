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
        // the plan gets no work area of its own, the Solver provides one shared buffer (see backend_*.cpp)
#if defined(BACKEND_CUDA)
    cufftType cuda_type = (xform == XFORM_R2C) ? CUFFT_D2Z : (xform == XFORM_C2R) ? CUFFT_Z2D
                                                                                  : CUFFT_Z2Z;
    backend_make_plan(&plan.fft, &plan.worksize, plan.dim, plan.idist, plan.odist, cuda_type, plan.howmany);
#elif defined(BACKEND_HIP)
    hipfftType hip_type = (xform == XFORM_R2C) ? HIPFFT_D2Z : (xform == XFORM_C2R) ? HIPFFT_Z2D
                                                                                   : HIPFFT_Z2Z;
    backend_make_plan(&plan.fft, &plan.worksize, plan.dim, plan.idist, plan.odist, hip_type, plan.howmany);
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
            // Pre-processing: fill the whole FFT input row (zero-padding included),
            // so the scratch buffer content is never read before being written
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % idist;
                const int j = idx / idist;

                const double w = (i == 0 || i == n - 1) ? 0.5 : 1.0;

                mem[j * idist + i] = (i < n) ? w * in[j * dist + stride * i] : 0.0;
            };
            parallel_for(pre_op, idist * howmany, "DCT Type I pre-processing", stream);

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
                const double re   = mem[2 * (j * odist + m)];
                const double im   = sign * mem[2 * (j * odist + m) + 1];

                double st, ct;  // sin/cos of pi*k/(2n)
                sincospi(k / (2.0 * n), &st, &ct);

                out[j * dist + stride * k] = 2.0 * (re * ct + im * st);
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
                    const double Xm  = in[j * dist + stride * m];
                    const double Xkc = in[j * dist + stride * kc];

                    double st, ct;  // sin/cos of pi*m/(2n)
                    sincospi(m / (2.0 * n), &st, &ct);

                    re = 0.5 * (ct * Xm + st * Xkc);
                    im = 0.5 * (st * Xm - ct * Xkc);
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
                    const int m = idx % n;
                    const int j = idx / n;

                    if (m < M) {
                        const double a = in[j * dist + stride * (2 * m)];
                        const double b = in[j * dist + stride * (2 * m + 1)];

                        double st, ct;  // sin/cos of pi*m/n
                        sincospi(m / (double)n, &st, &ct);

                        mem[2 * (j * idist + m)]     = a * ct - b * st;
                        mem[2 * (j * idist + m) + 1] = a * st + b * ct;
                    } else {
                        // zero-padding of the FFT input
                        mem[2 * (j * idist + m)]     = 0.0;
                        mem[2 * (j * idist + m) + 1] = 0.0;
                    }
                };
                parallel_for(pre_op, n * howmany, "DCT Type IV pre-processing", stream);

                flups_execute_c2c(plan, (flups_complex *)mem, (flups_complex *)mem, 1, stream);

                auto post_op = m_backend_lambda(const int idx) {
                    const int    k   = idx % n;
                    const int    j   = idx / n;
                    const double Cr  = mem[2 * (j * odist + k)];
                    const double Ci  = mem[2 * (j * odist + k) + 1];
                    const double C2r = mem[2 * (j * odist + (n - 1 - k))];
                    const double C2i = mem[2 * (j * odist + (n - 1 - k)) + 1];
                    const double Ca  = 0.5 * (Cr + C2r);
                    const double Sa  = 0.5 * (Ci - C2i);
                    const double Cb  = 0.5 * (Ci + C2i);
                    const double Sb  = 0.5 * (C2r - Cr);

                    double st, ct, s3t, c3t;  // sin/cos of theta = pi*(2k+1)/(4n) and 3*theta
                    sincospi((2 * k + 1) / (4.0 * n), &st, &ct);
                    sincospi(3.0 * (2 * k + 1) / (4.0 * n), &s3t, &c3t);

                    out[j * dist + stride * k] = 2.0 * (Ca * ct - Sa * st + Cb * c3t - Sb * s3t);
                };
                parallel_for(post_op, n * howmany, "DCT Type IV post-processing", stream);
            } else {
                // n odd fallback: complex pre-twiddle w[i]=x[i]*exp(i*pi*i/2n), a single
                // 2n-point complex FFT (inverse-direction kernel, sign=+1) of the
                // n-nonzero, n-zero-padded w, then a complex post-twiddle.
                auto pre_op = m_backend_lambda(const int idx) {
                    const int i = idx % dim;
                    const int j = idx / dim;

                    if (i < n) {
                        const double xi = in[j * dist + stride * i];

                        double st, ct;  // sin/cos of pi*i/(2n)
                        sincospi(i / (2.0 * n), &st, &ct);

                        mem[2 * (j * idist + i)]     = xi * ct;
                        mem[2 * (j * idist + i) + 1] = xi * st;
                    } else {
                        // zero-padding of the FFT input
                        mem[2 * (j * idist + i)]     = 0.0;
                        mem[2 * (j * idist + i) + 1] = 0.0;
                    }
                };
                parallel_for(pre_op, dim * howmany, "DCT Type IV pre-processing", stream);

                flups_execute_c2c(plan, (flups_complex *)mem, (flups_complex *)mem, 1, stream);

                auto post_op = m_backend_lambda(const int idx) {
                    const int    k  = idx % n;
                    const int    j  = idx / n;
                    const double re = mem[2 * (j * odist + k)];
                    const double im = mem[2 * (j * odist + k) + 1];

                    double st, ct;  // sin/cos of phi = pi*(2k+1)/(4n)
                    sincospi((2 * k + 1) / (4.0 * n), &st, &ct);

                    out[j * dist + stride * k] = 2.0 * (re * ct - im * st);
                };
                parallel_for(post_op, n * howmany, "DCT Type IV post-processing", stream);
            }

            break;
        }
        // DST Type I
        case FLUPS_DST_I: {
            // Pre-processing: fill the whole FFT input row (zero-padding included),
            // so the scratch buffer content is never read before being written
            auto pre_op = m_backend_lambda(const int idx) {
                const int i = idx % idist;
                const int j = idx / idist;

                mem[j * idist + i] = (1 <= i && i <= n) ? in[j * dist + stride * (i - 1)] : 0.0;
            };
            parallel_for(pre_op, idist * howmany, "DST Type I pre-processing", stream);

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
                const double re   = mem[2 * (j * odist + m)];
                const double im   = sign * mem[2 * (j * odist + m) + 1];

                double st, ct;  // sin/cos of pi*k/(2n)
                sincospi(k / (2.0 * n), &st, &ct);

                out[j * dist + stride * (n - 1 - k)] = 2.0 * (re * ct + im * st);
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
                    re = 0.5 * in[j * dist + stride * (n - 1)];
                    im = 0.0;
                } else {
                    const int    kc  = n - m;
                    const double Xm  = in[j * dist + stride * (n - 1 - m)];
                    const double Xkc = in[j * dist + stride * (n - 1 - kc)];

                    double st, ct;  // sin/cos of pi*m/(2n)
                    sincospi(m / (2.0 * n), &st, &ct);

                    re = 0.5 * (ct * Xm + st * Xkc);
                    im = 0.5 * (st * Xm - ct * Xkc);
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
                    const int m = idx % n;
                    const int j = idx / n;

                    if (m < M) {
                        const double a = in[j * dist + stride * (n - 1 - 2 * m)];
                        const double b = in[j * dist + stride * (n - 2 - 2 * m)];

                        double st, ct;  // sin/cos of pi*m/n
                        sincospi(m / (double)n, &st, &ct);

                        mem[2 * (j * idist + m)]     = a * ct - b * st;
                        mem[2 * (j * idist + m) + 1] = a * st + b * ct;
                    } else {
                        // zero-padding of the FFT input
                        mem[2 * (j * idist + m)]     = 0.0;
                        mem[2 * (j * idist + m) + 1] = 0.0;
                    }
                };
                parallel_for(pre_op, n * howmany, "DST Type IV pre-processing", stream);

                flups_execute_c2c(plan, (flups_complex *)mem, (flups_complex *)mem, 1, stream);

                auto post_op = m_backend_lambda(const int idx) {
                    const int    k    = idx % n;
                    const int    j    = idx / n;
                    const double Cr   = mem[2 * (j * odist + k)];
                    const double Ci   = mem[2 * (j * odist + k) + 1];
                    const double C2r  = mem[2 * (j * odist + (n - 1 - k))];
                    const double C2i  = mem[2 * (j * odist + (n - 1 - k)) + 1];
                    const double Ca   = 0.5 * (Cr + C2r);
                    const double Sa   = 0.5 * (Ci - C2i);
                    const double Cb   = 0.5 * (Ci + C2i);
                    const double Sb   = 0.5 * (C2r - Cr);
                    const double sign = (k % 2 == 0) ? 1.0 : -1.0;

                    double st, ct, s3t, c3t;  // sin/cos of theta = pi*(2k+1)/(4n) and 3*theta
                    sincospi((2 * k + 1) / (4.0 * n), &st, &ct);
                    sincospi(3.0 * (2 * k + 1) / (4.0 * n), &s3t, &c3t);

                    out[j * dist + stride * k] = sign * 2.0 * (Ca * ct - Sa * st + Cb * c3t - Sb * s3t);
                };
                parallel_for(post_op, n * howmany, "DST Type IV post-processing", stream);
            } else {
                auto pre_op = m_backend_lambda(const int idx) {
                    const int i = idx % dim;
                    const int j = idx / dim;

                    if (i < n) {
                        const double xi = in[j * dist + stride * (n - 1 - i)];

                        double st, ct;  // sin/cos of pi*i/(2n)
                        sincospi(i / (2.0 * n), &st, &ct);

                        mem[2 * (j * idist + i)]     = xi * ct;
                        mem[2 * (j * idist + i) + 1] = xi * st;
                    } else {
                        // zero-padding of the FFT input
                        mem[2 * (j * idist + i)]     = 0.0;
                        mem[2 * (j * idist + i) + 1] = 0.0;
                    }
                };
                parallel_for(pre_op, dim * howmany, "DST Type IV pre-processing", stream);

                flups_execute_c2c(plan, (flups_complex *)mem, (flups_complex *)mem, 1, stream);

                auto post_op = m_backend_lambda(const int idx) {
                    const int    k    = idx % n;
                    const int    j    = idx / n;
                    const double re   = mem[2 * (j * odist + k)];
                    const double im   = mem[2 * (j * odist + k) + 1];
                    const double sign = (k % 2 == 0) ? 1.0 : -1.0;

                    double st, ct;  // sin/cos of phi = pi*(2k+1)/(4n)
                    sincospi((2 * k + 1) / (4.0 * n), &st, &ct);

                    out[j * dist + stride * k] = sign * 2.0 * (re * ct - im * st);
                };
                parallel_for(post_op, n * howmany, "DST Type IV post-processing", stream);
            }

            break;
        }
        default:
            BACKEND_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }
}

#endif  // M_HAVE_GPU
