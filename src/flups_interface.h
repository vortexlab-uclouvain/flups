#ifndef FLUPS_INTERFACE_H_
#define FLUPS_INTERFACE_H_

/**
 * @brief to be used as "sign" for all of the FORWARD tranform
 *
 */
#define FLUPS_FORWARD -1  // equivalent to FFTW_FORWARD

/**
 * @brief to be used as "sign" for all of the BACKWARD tranform
 *
 */
#define FLUPS_BACKWARD 1  // equivalen to FFTW_BACKWARD

/**
 * @brief to be used as "sign" for all of the BACKWARD tranform
 *
 */
#define FLUPS_BACKWARD_DIFF 2

/**
 * @brief Memory alignment in bytes.
 *
 */
#define FLUPS_ALIGNMENT 16

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
/**
 * @name STRUCTURES AND DEFINITIONS
 * @{
 */
//==============================================================================

/**
 * @brief List of supported boundary conditions
 *
 * The boundary condition can be EVEN, ODD, PERiodic or UNBounded.
 */
enum BoundaryType {
    EVEN = 0, /**< EVEN boundary condition = zero flux  */
    ODD  = 1, /**< ODD boundary condition = zero value */
    PER  = 3, /**< PERiodic boundary conditions */
    UNB  = 4, /**< UNBounded boundary condition */
    NONE = 9  /**< No boundary condition = dimension not used */
};

/**
 * @brief The type of Green's function used for the Poisson solver
 *
 */
enum GreenType {
    CHAT_2 = 0, /**< @brief quadrature in zero, order 2, Chatelain et al. (2010) */
    LGF_2  = 1, /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2  = 2, /**< @brief regularized, order 2, Hejlesen et al. (2015)*/
    HEJ_4  = 3, /**< @brief regularized, order 4, Hejlesen et al. (2015)*/
    HEJ_6  = 4, /**< @brief regularized, order 6, Hejlesen et al. (2015)*/
    HEJ_8  = 5, /**< @brief regularized, order 8, Hejlesen et al. (2015)*/
    HEJ_10 = 6, /**< @brief regularized, order 10, Hejlesen et al. (2015)*/
    HEJ_0  = 7, /**< @brief Fourier cutoff, spectral-like, Hejlesen et al. (2019)*/
};

/**
 * @brief The type of possible solvers
 *
 * When solving for Biot-Savart, the Green's kernel \f$ G \f$ in Fourier space is adapted so that the Fourier
 * transform of the solution is obtained as
 *      \f[ \hat{\phi} = \hat{K} \times \hat{f} \f]
 * where \f$ \hat{K} \f$ is the spectral equivalent of the gradient of \f$ G \f$. Indeed,
 * the derivation is performed directly in Fourier space, according to the parameter @ref FLUPS_DiffType.
 */
enum SolverType {
    STD = 0, /**< @brief the standard poisson solver: \f$ \nabla^2(\phi) = (rhs) \f$ */
    ROT = 1  /**< @brief the Bio-Savart poisson solver: \f$ \nabla^2(\phi) = \nabla \times (rhs) \f$ */
};

/**
 * @brief The type of derivative to be used with @ref FLUPS_SolverType ROT.
 *
 */
enum DiffType {
    NOD = 0, /**< @brief Default parameter to be used with the STD type solve */
    SPE = 1, /**< @brief Spectral derivation, \f$ \hat{K} = i \, k \, \hat{G} \f$ */
    FD2 = 2, /**< @brief Spectral equivalent of 2nd order finite difference, \f$ \hat{K} = i \, \sin(kh) \, \hat{G} \f$ */
    FD4 = 4, /**< @brief Spectral equivalent of 4th order finite difference, \f$ \hat{K} = i \, ( 4/3 \sin(kh) - 1/6 \sin(2kh) ) \, \hat{G} \f$ */
    FD6 = 6  /**< @brief Spectral equivalent of 6th order finite difference, \f$ \hat{K} = i \, ( 3/2 \sin(kh) - 3/10 \sin(2kh) + 1/30 \sin(3kh) ) \, \hat{G} \f$ */
};

/**
 * @brief List of supported data center
 *
 * The data can be either Cell-centered or Node-centered
 */
enum CenterType {
    NODE_CENTER = 0, /**< NODE centered data  */
    CELL_CENTER = 1  /**< CELL centered data */
};

#endif
