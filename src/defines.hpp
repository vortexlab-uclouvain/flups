/**
 * @file defines.hpp
 * @copyright Copyright © Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <cassert>
#include <cmath>
#include <iostream>

#include <execinfo.h>
#include "fftw3.h"
#include "mpi.h"
#include "flups_interface.h" // get the main defines
#include "h3lpr/macros.hpp" // return the h3lpr macros
#include "h3lpr/profiler.hpp" // profiler
#include "h3lpr/ptr.hpp" // pointer allocation

//==============================================================================
//                      Compilation flags
//==============================================================================
/**
 * @brief FFTW planner flag driven by the NDEBUG flag
 *
 */
#ifndef FFTW_FLAG
#ifdef NDEBUG
#define FLUPS_FFTW_FLAG FFTW_PATIENT
#else
#define FLUPS_FFTW_FLAG FFTW_ESTIMATE
#endif
#else
#define FLUPS_FFTW_FLAG FFTW_FLAG
#endif

#ifndef COMM_DPREC
#define FLUPS_MPI_AGGRESSIVE 1
#else
#define FLUPS_MPI_AGGRESSIVE 0
#endif

/**
 * @brief enables the more evenly distributed balancing between ranks
 *
 * given N unknowns and P process, we try to evenly distribute the data
 * every rank has defacto B=N/P unknows as a baseline.
 * Then we are left with R = N%P unkowns to distribute among the P processes.
 * To do so instead of setting the R unknowns on the R first ranks we distribute them by groups.
 * We gather S (=stride) ranks together in a group and per group we add 1 unknow on the last rank of the group
 * The size of a group is given by S = P/R, which is also the stride between two groups
 * Exemple:
 *     - N = 32, P = 6: B = 5, R = 2 and therefore S = 3.
 *         So the rank distribution will be in two groups of 3 ranks:
 *              rank 0 -> 0 * 5 + 0 / 3 = 0
 *              rank 1 -> 1 * 5 + 1 / 3 = 5     (+5)
 *              rank 2 -> 2 * 5 + 2 / 3 = 10    (+5)
 *              rank 3 -> 3 * 5 + 3 / 3 = 16    (+6)
 *              rank 4 -> 4 * 5 + 4 / 3 = 21    (+5)
 *              rank 5 -> 5 * 5 + 5 / 3 = 26    (+5)
 *              rank 6 -> 6 * 5 + 6 / 3 = 32    (+6)
 *
 * To get the starting id from a rank we have:
 *       id = r * B + r/S
 *
 * To recover the rank from a global id (I) it's a bit longer.
 * We use S * B + 1 which is the number of unknowns inside one group
 * (1) get the group id:
 *      group_id = I /(S*B + 1)
 * (2) get the id within the group:
 *      local_group_id = I%(S*B + 1)
 * (3) get the rank within the group:
 *      local_group_id/B
 *
 * the rank is then:
 *      group_id * S + local_group_id/B
 */
#ifndef BALANCE_DPREC
#define FLUPS_NEW_BALANCE 1
#else
#define FLUPS_NEW_BALANCE 0
#endif

#ifdef HAVE_WISDOM
#define FLUPS_WISDOM_PATH HAVE_WISDOM
#endif

#ifdef HAVE_HDF5
#define FLUPS_HDF5 1
#else
#define FLUPS_HDF5 0
#endif

// register the current git commit for tracking purpose
#ifdef GIT_COMMIT
#define FLUPS_GIT_COMMIT GIT_COMMIT
#else
#define FLUPS_GIT_COMMIT "?"
#endif

#ifndef MPI_40
#define FLUPS_OLD_MPI 1
#else
#define FLUPS_OLD_MPI 0
#endif

#ifndef MPI_BATCH_SEND
#define FLUPS_MPI_BATCH_SEND 1
#else
#define FLUPS_MPI_BATCH_SEND MPI_BATCH_SEND
#endif

#ifndef MPI_MAX_NBSEND
#define FLUPS_MPI_MAX_NBSEND INT_MAX
#else
#define FLUPS_MPI_MAX_NBSEND MPI_MAX_NBSEND
#endif

#ifndef MPI_DEFAULT_ORDER
#define FLUPS_PRIORITYLIST 1
#else
#define FLUPS_PRIORITYLIST 0
#endif

#ifndef MPI_NO_ROLLING_RANK
#define FLUPS_ROLLING_RANK 1
#else
#define FLUPS_ROLLING_RANK 0
#endif

#ifndef MPI_NO_ALLOC
#define FLUPS_MPI_ALLOC 1
#else
#define FLUPS_MPI_ALLOC 0
#endif

//==============================================================================


#if (FLUPS_MPI_AGGRESSIVE)
#if (FLUPS_MPI_ALLOC)
using m_ptr_t = H3LPR::m_ptr<H3LPR::H3LPR_ALLOC_MPI, double*, FLUPS_ALIGNMENT>;
#else
using m_ptr_t = H3LPR::m_ptr<H3LPR::H3LPR_ALLOC_POSIX, double*, FLUPS_ALIGNMENT>;
#endif
#endif

/**
 * @brief allocates size bytes using the flups allocation (alignement)
 *
 * this funaction relies on the POSIX allocation as defined by H3LPR because
 * with a posix allocation the returned pointer is the one that needs to be freed
 *
 */
#define m_calloc(size)                                                             \
    ({                                                                             \
        H3LPR::m_ptr<H3LPR::H3LPR_ALLOC_POSIX, void *, FLUPS_ALIGNMENT> ptr(size); \
                                                                                   \
        ptr();                                                                     \
    })

/**
 * @brief free the memory allocated using m_calloc
 *
 */
#define m_free(data)     \
    {                    \
        std::free(data); \
    }



//=============================================================================
// Debug
//=============================================================================
#define FLUPS_print_data(topo, data)                                                                                                                                            \
    ({                                                                                                                                                                          \
        const int    ax0     = topo->axis();                                                                                                                                    \
        const int    ax1     = (ax0 + 1) % 3;                                                                                                                                   \
        const int    ax2     = (ax0 + 2) % 3;                                                                                                                                   \
        const int    lda     = topo->lda();                                                                                                                                     \
        const int    nf      = topo->nf();                                                                                                                                      \
        const int    nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};                                                                                                   \
        const size_t memdim  = topo->memdim();                                                                                                                                  \
        const size_t ondim   = topo->nloc(ax1) * topo->nloc(ax2);                                                                                                               \
        const size_t onmax   = topo->nloc(ax1) * topo->nloc(ax2);                                                                                                               \
        const size_t inmax   = topo->nloc(ax0);                                                                                                                                 \
        printf("ax0 = %d -- lda = %d -- nf = %d - nmem = %d %d %d -- end = %d %d %d \n", ax0, lda, nf, nmem[0], nmem[1], nmem[2], topo->nloc(0), topo->nloc(1), topo->nloc(2)); \
        for (int lia = 0; lia < lda; lia++) {                                                                                                                                   \
            printf("lia == %d \n", lia);                                                                                                                                        \
            if (nf == 1) {                                                                                                                                                      \
                for (int id = 0; id < onmax; id++) {                                                                                                                            \
                    const size_t   io     = id % ondim;                                                                                                                         \
                    opt_double_ptr argloc = data + lia * memdim + collapsedIndex(ax0, 0, io, nmem, nf);                                                                         \
                    if (id % topo->nloc(ax1) == 0) printf("\n");                                                                                                                \
                    for (size_t ii = 0; ii < inmax; ii++) {                                                                                                                     \
                        printf("%e \t ", argloc[ii]);                                                                                                                           \
                    }                                                                                                                                                           \
                    printf("\n");                                                                                                                                               \
                }                                                                                                                                                               \
            } else {                                                                                                                                                            \
                for (int id = 0; id < onmax; id++) {                                                                                                                            \
                    const size_t   io     = id % ondim;                                                                                                                         \
                    opt_double_ptr argloc = data + lia * memdim + collapsedIndex(ax0, 0, io, nmem, nf);                                                                         \
                    if (id % topo->nloc(ax1) == 0) printf("\n");                                                                                                                \
                    for (size_t ii = 0; ii < inmax; ii++) {                                                                                                                     \
                        printf("(%e, %e) \t", argloc[2 * ii], argloc[2 * ii + 1]);                                                                                              \
                    }                                                                                                                                                           \
                    printf("\n");                                                                                                                                               \
                }                                                                                                                                                               \
            }                                                                                                                                                                   \
        }                                                                                                                                                                       \
        fflush(stdout);                                                                                                                                                         \
    })
//=============================================================================
// WARNINGS
//=============================================================================
#define FLUPS_WARNING(format, ...)                   \
    ({                                             \
        m_log_def("FLUPS - Warning", format, ##__VA_ARGS__); \
    })


//=============================================================================
// LOG / FLUPS_INFO
//=============================================================================
#define FLUPS_INFO(format, ...)                   \
    ({                                              \
        m_verb_def("FLUPS ", format, ##__VA_ARGS__); \
    })

//=============================================================================
// VERBOSITY LEVELS
//=============================================================================
#if VERBOSE>=2
#define BEGIN_FUNC  m_begin_def("FLUPS")
#define END_FUNC m_end_def("FLUPS")
#else
#define BEGIN_FUNC \
    {              \
        ((void)0); \
    };
#define END_FUNC   \
    {              \
        ((void)0); \
    };
#endif

#if VERBOSE>=1
    #define FLUPS_INFO_1(...) FLUPS_INFO(__VA_ARGS__)
#else
    #define FLUPS_INFO_1(...) ((void)0);
#endif
#if VERBOSE>=2
    #define FLUPS_INFO_2(...) FLUPS_INFO(__VA_ARGS__)
#else
    #define FLUPS_INFO_2(...) ((void)0);
#endif
#if VERBOSE>=3
    #define FLUPS_INFO_3(...) FLUPS_INFO(__VA_ARGS__)
#else
    #define FLUPS_INFO_3(...) ((void)0);
#endif
#if VERBOSE>=4
    #define FLUPS_INFO_4(...) FLUPS_INFO(__VA_ARGS__)
#else
    #define FLUPS_INFO_4(...) ((void)0);
#endif

//=============================================================================
// ERRORS AND ASSERTS
//=============================================================================
#define FLUPS_CHECK(format, ...)                   \
    ({                                                \
        m_assert_def("FLUPS", format, ##__VA_ARGS__); \
    })


//=============================================================================
// CONSTANTS AND OTHERS
//=============================================================================

template <typename T>
static inline bool FLUPS_ISALIGNED(T a) {
    return ((uintptr_t)(const void*)a) % FLUPS_ALIGNMENT == 0;
}

template <typename T>
static inline int FLUPS_CMPT_ALIGNMENT(T a) {
    return ((uintptr_t)(const void*)a) % FLUPS_ALIGNMENT;
}

typedef int* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_int_ptr;
typedef double* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_double_ptr;
typedef fftw_complex* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_complex_ptr;

#define m_profStarti(prof, name, ...)                                    \
    ({                                                                   \
        H3LPR::Profiler *m_profStarti_prof_ = (H3LPR::Profiler *)(prof); \
        char m_profStarti_name_[1024];                                   \
        std::sprintf(m_profStarti_name_, name, ##__VA_ARGS__);           \
        m_profStart(m_profStarti_prof_, m_profStarti_name_);             \
    })

#define m_profStopi(prof, name, ...)                                    \
    ({                                                                  \
        H3LPR::Profiler *m_profStopi_prof_ = (H3LPR::Profiler *)(prof); \
        char m_profStopi_name_[1024];                                   \
        sprintf(m_profStopi_name_, name, ##__VA_ARGS__);                \
        m_profStop(m_profStopi_prof_, m_profStopi_name_);               \
    })

//=============================================================================
// MEMORY ALLOCATION AND FREE
//=============================================================================

#if defined(__INTEL_COMPILER)
    #define FLUPS_ASSUME_ALIGNED(a,b) __assume_aligned(a,b)
#elif defined(__GNUC__)
    #define FLUPS_ASSUME_ALIGNED(a,b) __builtin_assume_aligned(a,b)
#endif

// typedef enum FLUPS_BoundaryType BoundaryType;
// typedef enum FLUPS_GreenType    GreenType;

static const double c_1opi     = 1.0 / (1.0 * M_PI);
static const double c_1o2pi    = 1.0 / (2.0 * M_PI);
static const double c_1o2pi2   = 1.0 / (2.0 * M_PI*M_PI);
static const double c_1o4pi    = 1.0 / (4.0 * M_PI);
static const double c_1osqrtpi = 1.0 / sqrt(M_PI);
static const double c_sqrtpi   = sqrt(M_PI);
static const double c_1o2      = 1.0 / 2.0;
static const double c_1o4      = 1.0 / 4.0;
static const double c_7o4      = 7. / 4.;
static const double c_19o8     = 19. / 8;
static const double c_2o3      = 2. / 3;
static const double c_1o12     = 1. / 12.;
static const double c_1o24     = 1. / 24;
static const double c_5o16     = 5. / 16.;
static const double c_1o16     = 1. / 16.;
static const double c_187o64   = 187. / 64.;
static const double c_233o192  = 233. / 192.;
static const double c_29o192   = 29. / 192.;
static const double c_13o192   = 13. / 192.;
static const double c_1o192    = 1. / 192.;
static const double c_11o12    = 11. / 12.;
static const double c_7o24     = 7. / 24.;
static const double c_25o24    = 25. / 24.;
static const double c_1o48     = 1. / 48.;
static const double c_23o48    = 23. / 48.;
static const double c_1o384    = 1. / 384.;
static const double c_17o384   = 17. / 384.;
static const double c_23o384   = 23. / 384.;
static const double c_11o32    = 11. / 32.;
static const double c_1o96     = 1. / 96.;
static const double c_35o128   = 35. / 128.;
static const double c_47o128   = 47. / 128.;
static const double c_93o256   = 93. / 256.;
static const double c_47o256   = 47. / 256.;
static const double c_1o256    = 1. / 256.;
static const double c_73o768   = 73. / 768.;
static const double c_11o768   = 11. / 768.;
static const double c_23o768   = 23. / 768.;
static const double c_1o768    = 1. / 768.;
static const double c_4o3      = 4. / 3.;
static const double c_1o6      = 1. / 6.;
static const double c_3o2      = 3. / 2.;
static const double c_3o10     = 3. / 10.;
static const double c_1o30     = 1. / 30.;

static const double c_1osqrt2 = 1.0 / M_SQRT2;

static const double c_2pi = 2.0 * M_PI;

#endif