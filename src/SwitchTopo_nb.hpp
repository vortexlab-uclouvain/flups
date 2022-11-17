/**
 * @file SwitchTopo_nb.hpp
 * @copyright Copyright © Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef SWITCHTOPO_NB_HPP
#define SWITCHTOPO_NB_HPP

#include <cstring>
#include "defines.hpp"
#include "hdf5_io.hpp"
#include "mpi.h"
#include "Topology.hpp"
#include "omp.h"
#include "SwitchTopo.hpp"

/**
 * @brief Switch between to different topologies using non-blocking communications
 * 
 * Reorganize the memory between 2 different topologies, also accounting for a
 * "principal axis" which is aligned with the fast rotating index.
 * 
 * Communications are handled by packing data in blocks (i.e. chuncks). These 
 * smaller pieces are sent and recieved asynchronesouly, in a hope for reducing
 * the time required for communication.
 * 
 * The a memory shift can be specified in the switch between the topologies, when
 * there is a need for skipping some data at the left or right side of a given direction.
 * 
 */
class SwitchTopo_nb : public SwitchTopo
{
protected:
    int selfBlockN_ = 0;
    int *iselfBlockID_ = NULL; /**<@brief The list of the block iD that stays on the current rank in the input topology (used while output to input)  */
    int *oselfBlockID_ = NULL; /**<@brief The list of the block iD that stays on the current rank in the output topoloy (used while input to ouput) */

    int *i2o_destTag_ = NULL; /**<@brief The destination rank in the output topo of each block */
    int *o2i_destTag_ = NULL; /**<@brief The destination rank in the output topo of each block */

    MPI_Request *i2o_sendRequest_ = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *i2o_recvRequest_ = NULL; /**<@brief The MPI Request generated on the recv */
    MPI_Request *o2i_sendRequest_ = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *o2i_recvRequest_ = NULL; /**<@brief The MPI Request generated on the recv */

    void init_blockInfo_(const Topology *topo_in, const Topology *topo_out);
    void free_blockInfo_();

   public:
    SwitchTopo_nb(const Topology *topo_input, const Topology *topo_output, const int shift[3],H3LPR::Profiler* prof);
    ~SwitchTopo_nb();

    void setup_buffers(opt_double_ptr sendBuf_,opt_double_ptr recvBuf_);
    void execute(double* v, const int sign) const;
    void setup() ;
    void disp() const;
};

void SwitchTopo_test();

#endif