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
    CHAT_2 = 0,     /**< @brief quadrature in zero, order 2, Chatelain et al. (2010) */
    LGF_2  = 1,     /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2  = 2,     /**< @brief regularized, order 2, Hejlesen et al. (2015)*/
    HEJ_4  = 3,     /**< @brief regularized, order 4, Hejlesen et al. (2015)*/
    HEJ_6  = 4,     /**< @brief regularized, order 6, Hejlesen et al. (2015)*/
    HEJ_8  = 5,     /**< @brief regularized, order 8, Hejlesen et al. (2015)*/
    HEJ_10 = 6,     /**< @brief regularized, order 10, Hejlesen et al. (2015)*/
    HEJ_0  = 7,     /**< @brief Fourier cutoff, spectral-like, Hejlesen et al. (2019)*/
    LGF_4  = 8,     /**< @brief Lattice Green's function, order 4 */
    LGF_6  = 9,     /**< @brief Lattice Green's function, order 6 */
    LGF_8  = 10,    /**< @brief Lattice Green's function, order 8 */
    MEHR_4L = 11,   /**< @brief Spectral equivalent of the left-hand side of the Mehrstellen HOC stencil, order 4 */
    MEHR_6L = 12,   /**< @brief Spectral equivalent of the left-hand side of the Mehrstellen HOC stencil, order 6 */
    MEHR_4F = 13,   /**< @brief Spectral equivalent of the full Mehrstellen HOC stencil, order 4 */
    MEHR_6F = 14,   /**< @brief Spectral equivalent of the full Mehrstellen HOC stencil, order 6 */
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
