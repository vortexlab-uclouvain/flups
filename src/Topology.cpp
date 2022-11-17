/**
 * @file Topology.cpp
 * @copyright Copyright © Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#include "Topology.hpp"
#include "hdf5_io.hpp"


/**
 * @brief Construct a new Topology
 * 
 * @param nf the number of field (eg: real = 1, complex = 2)
 * @param lda leading dimension of array=the number of components (eg scalar=1, vector=3)
 * @param axis the dimension of the fastest rotating index (eg: 0, 1 or 2)
 * @param nglob the global size per dim
 * @param nproc the number of proc per dim
 * @param isComplex indicate if the Topology uses complex indexing or not
 * @param axproc gives the order of the rank decomposition (eg. (0,2,1) to start decomposing in X then Z then Y). If NULL is passed, use by default (0,1,2).
 * @param alignment the number of bytes on which we want the topology to be aligned along the #axis only
 * @param comm the communicator associated to the topology.
 * 
 * If the MPI comm is associated with a MPI_CART topology, axproc is ignored and we use the MPI routines to determine the 3D rank from the global rank (and vice versa).
 * 
 */
Topology::Topology(const int axis, const int lda, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment, MPI_Comm comm):alignment_(alignment) {
    BEGIN_FUNC;

    comm_ = comm;
    lda_  = lda;

    int comm_size, rank;
    MPI_Comm_size(comm_,&comm_size);
    MPI_Comm_rank(comm_,&rank);

    FLUPS_CHECK(nproc[0]*nproc[1]*nproc[2] == comm_size,"the total number of procs (=%d) have to be = to the comm size (=%d)",nproc[0]*nproc[1]*nproc[2], comm_size);

    //-------------------------------------------------------------------------
    /** - get memory axis and complex information  */
    //-------------------------------------------------------------------------
    axis_ = axis;
    if (!isComplex) {
        nf_ = 1;
    } else {
        nf_ = 2;
    }

    //-------------------------------------------------------------------------
    /** - get proc information  */
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        // total dimension
        nglob_[id] = nglob[id];
        // number of proc in each direction
        nproc_[id] = nproc[id];
        // store the proc axis, used to split the rank
        axproc_[id] = (axproc == NULL) ? id : axproc[id];
    }

    //-------------------------------------------------------------------------
    /** - split the rank and get rankd  */
    //-------------------------------------------------------------------------
    ranksplit(rank, axproc_, nproc_, comm_, rankd_);

    //-------------------------------------------------------------------------
    /** - Get the new sizes  */
    //-------------------------------------------------------------------------
    cmpt_sizes();

    FLUPS_INFO("New topo created: nf = %d, axis = %d, local sizes = %d %d %d vs mem size = %d %d %d -- global sizes = %d %d %d",nf_,axis_,nloc_[0],nloc_[1],nloc_[2],nmem_[0],nmem_[1],nmem_[2],nglob_[0],nglob_[1],nglob_[2]);
    END_FUNC;
}

/**
 * @brief compute the nloc and nmem sizes using rankd_, nglob_, nproc_, nloc_
 * 
 * This function padds the size of the domain if needed
 * 
 */
void Topology::cmpt_sizes() {
    BEGIN_FUNC;
    for (int id = 0; id < 3; id++) {
        // we get the max between the nglob and
        nloc_[id] = cmpt_nbyproc(id);
        nmem_[id] = nloc_[id];
        // if we are in the axis and the last proc, we pad to ensure that every pencil is ok with alignment
        // if (id == axis_ && rankd_[id] == (nproc_[id] - 1)) {
        if (id == axis_) {
            // compute by how many we are not aligned: the global size in double = nglob * nf
            const int modulo = (nloc_[id] * nf_ * sizeof(double)) % alignment_;
            // compute the number of points to add (in double indexing)
            const int delta = (alignment_ - modulo) / sizeof(double);
            nmem_[id] += (modulo == 0) ? 0 : delta / nf_;
        }
    }
    END_FUNC;
}

/**
 * @brief Set a new communicator for the topology
 * 
 * @param comm 
 */
void Topology::change_comm(MPI_Comm comm) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - Store the new communicator  */
    //-------------------------------------------------------------------------
    comm_ = comm;
    //-------------------------------------------------------------------------
    /** - Recompute the new rank and rankd info  */
    //-------------------------------------------------------------------------
    int newrank;
    MPI_Comm_rank(comm, &newrank);
    // split the rank to get the new rankd
    ranksplit(newrank, axproc_, nproc_, comm_, rankd_);
    //-------------------------------------------------------------------------
    /** - Get the sizes  */
    //-------------------------------------------------------------------------
    cmpt_sizes();

    END_FUNC;
}

/**
 * @brief Destroy the Topology
 * 
 */
Topology::~Topology() {}

/**
 * @brief compute the starting and ending ids for the current topo in order to be inside the other topology's bounds
 *
 * @param shift the shift between the 2 topos: current topo in (0,0,0) = other topo in (shift)
 * @param other the other topology
 * @param start the start index to use in the current topo (in the local indexing!)
 * @param end the end index to use in the current topo (in the local indexing!)
 */
void Topology::cmpt_intersect_id(const int shift[3], const Topology* other, int start[3], int end[3]) const {
    BEGIN_FUNC;
    FLUPS_CHECK(this->isComplex() == other->isComplex(), "The two topo have to be both complex or real");
    //--------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        // get the maximum index in the other topology, minimum is 0 by definition
        const int onglob = other->nglob(id);

        // we try to get the first point
        const int start_global_id = shift[id] + this->cmpt_start_id(id);
        const int end_global_id   = shift[id] + this->cmpt_start_id(id) + this->nloc(id);

        // enforces the bounds
        // if the start_global_id is negative, we have offset the start index
        start[id] = -m_min(start_global_id, 0);
        // if the end_global_id is bigger than onglob we have to offset the end
        end[id] = this->nloc(id) - m_max(end_global_id - onglob, 0);
    }
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief display topology most important info
 * 
 */
void Topology::disp() const {
    BEGIN_FUNC;

    int comm_size, rank;
    MPI_Comm_size(comm_,&comm_size);
    MPI_Comm_rank(comm_,&rank);

    FLUPS_INFO("------------------------------------------");
    FLUPS_INFO("## Topology created on proc %d/%d", rank, comm_size);
    FLUPS_INFO(" - axis = %d", axis_);
    FLUPS_INFO(" - lda = %d", lda_);
    FLUPS_INFO(" - nglob = %d %d %d", nglob_[0], nglob_[1], nglob_[2]);
    FLUPS_INFO(" - nloc = %d %d %d", nloc_[0], nloc_[1], nloc_[2]);
    FLUPS_INFO(" - nmem = %d %d %d", nmem_[0], nmem_[1], nmem_[2]);
    FLUPS_INFO(" - nproc = %d %d %d", nproc_[0], nproc_[1], nproc_[2]);
    FLUPS_INFO(" - rankd = %d %d %d", rankd_[0], rankd_[1], rankd_[2]);
    FLUPS_INFO(" - axproc = %d %d %d", axproc_[0], axproc_[1], axproc_[2]);
    FLUPS_INFO(" - isComplex = %d", nf_ == 2);
    FLUPS_INFO("------------------------------------------");
}

void Topology::disp_rank() {
    BEGIN_FUNC;
#ifdef DUMP_DBG
    // we only focus on the real size = local size
    double* rankdata = (double*)m_calloc(sizeof(double) * this->locsize() * 2);
    int     rank, rank_new;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_rank(comm_, &rank_new);

    for (int i = 0; i < this->locsize(); i++) {
        rankdata[2 * i]     = rank + rank_new / 100.;
        rankdata[2 * i + 1] = rankd_[0] + rankd_[1] / 10. + rankd_[2] / 100.;
    }

    int  rlen;
    char commname[MPI_MAX_OBJECT_NAME];
    MPI_Comm_get_name(comm_, commname, &rlen);
    std::string cn(commname, rlen);
    std::string name = "rank_topo_axis" + std::to_string(this->axis()) + "procs_" + std::to_string(this->nproc(0)) + std::to_string(this->nproc(1)) + std::to_string(this->nproc(2)) + "_" + cn;
    if (this->isComplex()) {
        hdf5_dump(this, name, rankdata);
    } else {
        this->switch2complex();
        hdf5_dump(this, name, rankdata);
        this->switch2real();
    }
    m_free(rankdata);
#endif
    END_FUNC;
}