/**
 * @file SwitchTopo.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-25
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#ifndef REODER_MPI_HPP
#define REODER_MPI_HPP

#include <cstring>
#include "defines.hpp"
#include "hdf5_io.hpp"
#include "mpi.h"
#include "topology.hpp"

/**
 * @brief Do the switch between to different topologies
 * 
 * The switch between two topologies is driven by the shift between the topologies.
 * 
 */
class SwitchTopo {
   protected:
    const Topology *_topo_in;  /**<@brief input topology  */
    const Topology *_topo_out; /**<@brief  output topology */

    int _ishift[3]; /**<@brief the shift index for #_topo_in to match coordinate (0,0,0) in #__topo_out   */
    int _oshift[3]; /**<@brief the shift index for #_topo_out to match coordinate (0,0,0) in #_topo_in*/

    int _istart[3]; /**<@brief the starting index for #_topo_in to be inside #_topo_out  */
    int _iend[3];   /**<@brief the ending index for #_topo_in to be inside #_topo_out  */
    int _ostart[3]; /**<@brief the starting index for #_topo_out to be inside #_topo_in  */
    int _oend[3];   /**<@brief the ending index for #_topo_out to be inside #_topo_in  */

    int *_nsend          = NULL; /**<@brief number of unknowns to send to each proc  */
    int *_nrecv          = NULL; /**<@brief number of unknowns to receiv from each proc */
    int *_ssend          = NULL; /**<@brief start index in the buffer memory for data to send to each proc */
    int *_srecv          = NULL; /**<@brief start index in the buffer memory for data to receive from each proc */
    int *_count          = NULL; /**<@brief counter to know where to put the current readen data  */
    int *_naxis_i2o_send = NULL; /**<@brief naxis when sending from the _topo_in to _topo_out  */
    int *_naxis_o2i_send = NULL; /**<@brief naxis when sending from the _topo_out to _topo_in */
    int *_naxis_i2o_recv = NULL; /**<@brief naxis when receiving from the _topo_in to _topo_out  */
    int *_naxis_o2i_recv = NULL; /**<@brief naxis when receiving from the _topo_out to _topo_in */

    double *_bufsend = NULL; /**<@brief buffer to send data */
    double *_bufrecv = NULL; /**<@brief buffer to receive data */

   public:
    SwitchTopo(const Topology *topo_input, const Topology *topo_output, const int shift[3]);
    ~SwitchTopo();

    void execute(opt_double_ptr v, const int sign);

    void disp();
    void test();
};

#endif