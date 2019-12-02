/**
 * @file hdf5_io.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
 * 
 */

#include "hdf5_io.hpp"

void hdf5_dump(const Topology *topo, const string filename, const double *data) {
    BEGIN_FUNC;
    xmf_write(topo, filename, "data");
    hdf5_write(topo, filename, "data", data);
    END_FUNC;
}

/**
 * @brief writes the hdf5 file referenced in an xmf file
 * 
 * If the topology is real it creates one dataset `attribute`.
 * If the topology is complex it creates two datasets `attribute_real` and `attribute_imag`
 * 
 * 
 * @param topo topology of data being exported
 * @param filename the filename for export without its extension
 * @param attribute the name of the attribute to be added
 * @param data the array containing the data associated to the #topo
 * 
 * ---------------------------------------
 * Inspired from https://portal.hdfgroup.org/display/HDF5/Writing+by+Chunk+in+PHDF5
 * 
 * **General comments:**
 * A dataset is build on a dataspace (describing the datalayout) and a datatype (describing one element)
 * 
 * @warning
 * - the fastest rotating index is located at the last positions in the hdf5 C convention: `(axis+2 , axis+1 , axis)`
 * - the stride in the hyperslabs functions is the stride of 2 successive blocks!!
 * 
 * **Compression:**
 * The compression used is GNU gzip. The compression is characterize by the level of compression:
 * 
 * - `0`: No compression
 * - `1`: Best compression speed, least compression ratio
 * - `2` to `8`: compression speed degrades, compression ratio improves
 * - `9`: Slowest compression speed, best compression ratio
 * 
 * 
 * @todo for the moment no compression is set since it fails when changing the type from double to float
 * 
 * We do the following steps:
 * 
 */
void hdf5_write(const Topology *topo, const string filename, const string attribute, const double *data) {
    BEGIN_FUNC;

    int mpi_size, mpi_rank;
    MPI_Comm comm = topo->get_comm();
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    hid_t  file_id;                         // file id
    hid_t  filespace_real, filespace_imag;  //dataspaces
    hid_t  fileset_real, fileset_imag;      // datasets
    hid_t  memspace;
    hid_t  plist_id; /* property list identifier */
    herr_t status;   // error code

    string extFilename = "data/" + filename + ".h5";

    // compression level
    // int complevel = 0;

    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    //-------------------------------------------------------------------------
    /** - Create a new file collectively  */
    //-------------------------------------------------------------------------
    // setup the property list for file access (property list = option list)
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    // do the magic stuff
    MPI_Info FILE_INFO_TEMPLATE;
    MPI_Info_create(&FILE_INFO_TEMPLATE);
    H5Pset_sieve_buf_size(plist_id, 262144);
    H5Pset_alignment(plist_id, 524288, 262144);
    MPI_Info_set(FILE_INFO_TEMPLATE, "access_style", "write_once");
    MPI_Info_set(FILE_INFO_TEMPLATE, "collective_buffering", "true");
    MPI_Info_set(FILE_INFO_TEMPLATE, "cb_block_size", "1048576");
    MPI_Info_set(FILE_INFO_TEMPLATE, "cb_buffer_size", "4194304");
    // do some magic
    // H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);
    H5Pset_fapl_mpio(plist_id, comm, FILE_INFO_TEMPLATE);
    MPI_Info_free(&FILE_INFO_TEMPLATE);
    // create the file ID
    file_id = H5Fcreate(extFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    if (file_id < 0) FLUPS_ERROR("Failed to open the file.", LOCATION);
    // close the property list
    H5Pclose(plist_id);

    //-------------------------------------------------------------------------
    /** - Create the file dataspace and dataset  */
    //-------------------------------------------------------------------------
    // the file information is given by the global size

    hsize_t* field_dims  = (hsize_t*) malloc( ((int) (topo->lda()!=1) + 3)*sizeof(hsize_t) );
    if(topo->lda()==1){

    field_dims[0] = (hsize_t)topo->nglob(ax2);  // field global dimensions
    field_dims[1] = (hsize_t)topo->nglob(ax1);  // field global dimensions
    field_dims[2] = (hsize_t)topo->nglob(ax0);  // field global dimensions
    }else{
    field_dims[0] = topo->lda();
    field_dims[1] = (hsize_t)topo->nglob(ax2);  // field global dimensions
    field_dims[2] = (hsize_t)topo->nglob(ax1);  // field global dimensions
    field_dims[3] = (hsize_t)topo->nglob(ax0);  // field global dimensions
    }

    // setup the property list = option list
    plist_id           = H5Pcreate(H5P_DATASET_CREATE);
    if(topo->lda()==1){
    hsize_t chk_dim[3] = {8, 8, 8};
    H5Pset_chunk(plist_id, 3, chk_dim);
    }else{
     hsize_t chk_dim[4] = {8, 8, 8, 8};
    H5Pset_chunk(plist_id, 4, chk_dim);
    }

    // create dataset and dataspace
    //
    printf("LDA? %d\n",(int) (topo->lda()!=1) + 3);
    if (!topo->isComplex()) {
        filespace_real = H5Screate_simple((int) (topo->lda()!=1) + 3, field_dims, NULL);
        fileset_real   = H5Dcreate(file_id, attribute.c_str(), H5T_NATIVE_FLOAT, filespace_real, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Sclose(filespace_real);
    } else {
        string realname = attribute + "_real";
        filespace_real  = H5Screate_simple((int) (topo->lda()!=1) + 3, field_dims, NULL);
        fileset_real    = H5Dcreate(file_id, realname.c_str(), H5T_NATIVE_FLOAT, filespace_real, H5P_DEFAULT, plist_id, H5P_DEFAULT);

        string imagname = attribute + "_imag";
        filespace_imag  = H5Screate_simple((int) (topo->lda()!=1) + 3, field_dims, NULL);
        fileset_imag    = H5Dcreate(file_id, imagname.c_str(), H5T_NATIVE_FLOAT, filespace_imag, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    }
    // close property list
    H5Pclose(plist_id);

    free(field_dims);

    //-------------------------------------------------------------------------
    /** - select the hyperslab inside the file dataset (=writting location)  */
    //-------------------------------------------------------------------------
    // get the offset from topo
    int topo_offset[3];
    topo->get_istart_glob(topo_offset);

    hsize_t* count  = (hsize_t*) malloc( ((int) (topo->lda()!=1) + 3)*sizeof(hsize_t) );
    hsize_t* stride = (hsize_t*) malloc( ((int) (topo->lda()!=1) + 3)*sizeof(hsize_t) );
    hsize_t* block  = (hsize_t*) malloc( ((int) (topo->lda()!=1) + 3)*sizeof(hsize_t) );
    hsize_t* offset = (hsize_t*) malloc( ((int) (topo->lda()!=1) + 3)*sizeof(hsize_t) );


    for (int i = 0; i<3; i++){
        count[i]  = 1;                                                                          // how many blocks to write
        stride[i] = 1;                                                                          // distance between 2 blocks
    }


        // compute some memory quantities
    if (topo->lda() == 1) {
        block[0]  = (hsize_t)topo->nloc(ax2);   // the block size = the local size
        block[1]  = (hsize_t)topo->nloc(ax1);   // the block size = the local size
        block[2]  = (hsize_t)topo->nloc(ax0);   // the block size = the local size
        offset[0] = (hsize_t)topo_offset[ax2];  // offset in the file
        offset[1] = (hsize_t)topo_offset[ax1];  // offset in the file
        offset[2] = (hsize_t)topo_offset[ax0];  // offset in the file
    } else {
        count[3]  = 1;  // how many blocks to write
        stride[3] = 1;  // distance between 2 blocks
        block[0]  = topo->lda();
        block[1]  = (hsize_t)topo->nloc(ax2);  // the block size = the local size
        block[2]  = (hsize_t)topo->nloc(ax1);  // the block size = the local size
        block[3]  = (hsize_t)topo->nloc(ax0);  // the block size = the local size
        offset[0] = 0;
        offset[1] = (hsize_t)topo_offset[ax2];  // offset in the file
        offset[2] = (hsize_t)topo_offset[ax1];  // offset in the file
        offset[3] = (hsize_t)topo_offset[ax0];  // offset in the file
    }

    // get the hyperslab within the dataset
    if (!topo->isComplex()) {
        filespace_real = H5Dget_space(fileset_real);
        status         = H5Sselect_hyperslab(filespace_real, H5S_SELECT_SET, offset, stride, count, block);
        FLUPS_CHECK(status >= 0, "Failed to select hyperslab in dataset.", LOCATION);
    } else {
        filespace_real = H5Dget_space(fileset_real);
        status         = H5Sselect_hyperslab(filespace_real, H5S_SELECT_SET, offset, stride, count, block);
        FLUPS_CHECK(status >= 0, "Failed to select real hyperslab in dataset.", LOCATION);
        filespace_imag = H5Dget_space(fileset_imag);
        status         = H5Sselect_hyperslab(filespace_imag, H5S_SELECT_SET, offset, stride, count, block);
        FLUPS_CHECK(status >= 0, "Failed to select complex hyperslab in dataset.", LOCATION);
    }

    free(count);
    free(stride);
    free(block);
    free(offset);

    //-------------------------------------------------------------------------
    /** - do the writting  */
    //-------------------------------------------------------------------------
    // dataspace = data inside the memory that has a full size of nmem
    if(topo->lda()==1){
    hsize_t memsize[3] = {(hsize_t)topo->nmem(ax2), (hsize_t)topo->nmem(ax1), (hsize_t)(topo->nmem(ax0) * topo->nf())};
    memspace           = H5Screate_simple(3, memsize, NULL);
    } else {
    hsize_t memsize[4] = {topo->lda(), (hsize_t)topo->nmem(ax2), (hsize_t)topo->nmem(ax1), (hsize_t)(topo->nmem(ax0) * topo->nf())};
    memspace           = H5Screate_simple(4, memsize, NULL);
 }

    // set the property list
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    hsize_t* memcount = (hsize_t*) malloc( ((int) (topo->lda()!=1) + 3)*sizeof(hsize_t) );
    hsize_t* memblock = (hsize_t*) malloc( ((int) (topo->lda()!=1) + 3)*sizeof(hsize_t) );

    // get the data counts that we are going to write -> restricted to nloc instead of nmem
    if(topo->lda()==1){
	    memblock[0] = 1;
	    memblock[1] = 1;
	    memblock[2] = 1;
            memcount[0] = (hsize_t)topo->nloc(ax2);
            memcount[1] = (hsize_t)topo->nloc(ax1);
            memcount[2] = (hsize_t)topo->nloc(ax0);
    }else {
            memblock[0] = 1;
	    memblock[1] = 1;
	    memblock[2] = 1;
	    memblock[3] = 1;
	    memcount[0] = topo->lda();
	    memcount[1] = (hsize_t)topo->nloc(ax2);
	    memcount[2] = (hsize_t)topo->nloc(ax1);
	    memcount[3] = (hsize_t)topo->nloc(ax0);
    }

    if (!topo->isComplex()) {
        // get the hyperslab within the dataset
        if (topo->lda() == 1) {
            hsize_t memoffset[3] = {0, 0, 0};  // offset in memory
            hsize_t memstride[3] = {1, 1, 1};
            status               = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
            FLUPS_CHECK(status >= 0, "Failed to select hyperslab in memmory.", LOCATION);
            status = H5Dwrite(fileset_real, H5T_NATIVE_DOUBLE, memspace, filespace_real, plist_id, data);
            FLUPS_CHECK(status >= 0, "Failed to write hyperslab to file.", LOCATION);
        } else {
            hsize_t memoffset[4] = {0, 0, 0, 0};  // offset in memory
            hsize_t memstride[4] = {1, 1, 1, 1};
            status               = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
            FLUPS_CHECK(status >= 0, "Failed to select hyperslab in memmory.", LOCATION);
            status = H5Dwrite(fileset_real, H5T_NATIVE_DOUBLE, memspace, filespace_real, plist_id, data);
            FLUPS_CHECK(status >= 0, "Failed to write hyperslab to file.", LOCATION);
        }
    }

    if (topo->isComplex()) {
        if(topo->lda()==1){
            // stride is 2 for complex numbers
            hsize_t memstride[3] = { 1, 1, 2};
            // real part
            hsize_t memoffset[3] = { 0, 0, 0};
            status               = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
            FLUPS_CHECK(status >= 0, "Failed to select real hyperslab in memmory.", LOCATION);
            status = H5Dwrite(fileset_real, H5T_NATIVE_DOUBLE, memspace, filespace_real, plist_id, data);
            FLUPS_CHECK(status >= 0, "Failed to write real part hyperslab to file.", LOCATION);

            memoffset[2] = 1;  // set an offset on the fastest rotating index

                status       = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
            FLUPS_CHECK(status >= 0, "Failed to select imag hyperslab in memmory.", LOCATION);
            status = H5Dwrite(fileset_imag, H5T_NATIVE_DOUBLE, memspace, filespace_imag, plist_id, data);
            FLUPS_CHECK(status >= 0, "Failed to write imaginary part hyperslab to file.", LOCATION);

    
        } else{
            // stride is 2 for complex numbers
            hsize_t memstride[4] = {1, 1, 1, 2};
            // real part
            hsize_t memoffset[4] = {0, 0, 0, 0};
            status               = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
            FLUPS_CHECK(status >= 0, "Failed to select real hyperslab in memmory.", LOCATION);
            status = H5Dwrite(fileset_real, H5T_NATIVE_DOUBLE, memspace, filespace_real, plist_id, data);
            FLUPS_CHECK(status >= 0, "Failed to write real part hyperslab to file.", LOCATION);
        
            memoffset[3] = 1;  // set an offset on the fastest rotating index

            // imaginary part

            status       = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
            FLUPS_CHECK(status >= 0, "Failed to select imag hyperslab in memmory.", LOCATION);
            status = H5Dwrite(fileset_imag, H5T_NATIVE_DOUBLE, memspace, filespace_imag, plist_id, data);
            FLUPS_CHECK(status >= 0, "Failed to write imaginary part hyperslab to file.", LOCATION);
        }
    }

free(memcount);
free(memblock);

    //-------------------------------------------------------------------------
    /** - close everything */
    //-------------------------------------------------------------------------
    H5Dclose(fileset_real);
    H5Sclose(filespace_real);
    H5Sclose(memspace);
    if (topo->isComplex()) {
        H5Dclose(fileset_imag);
        H5Sclose(filespace_imag);
    }

    H5Pclose(plist_id);
    H5Fclose(file_id);

    END_FUNC;
    return;
}

/**
 * @brief writes a xmf file, readable by a XDMF viewer
 * 
 * @param topo topology of data being exported
 * @param filename the filename for export without its extension
 * @param attribute the name of the attribute to be added
 * 
 * -----------------------------
 * inspired from http://visitusers.org/index.php?title=Using_XDMF_to_read_HDF5
 * see also http://www.xdmf.org/index.php/XDMF_Model_and_Format
 */
void xmf_write(const Topology *topo, const string filename, const string attribute) {
    BEGIN_FUNC;

    int rank;
    MPI_Comm comm = topo->get_comm();
    MPI_Comm_rank(comm, &rank);

    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    FILE * xmf         = 0;

    string folder = "./data";
    string extFilename = folder + "/" + filename + ".xmf";

    if (rank == 0) {
        struct stat st = {0};
        if (stat(folder.c_str(), &st) == -1) {
                mkdir(folder.c_str(), 0770); //create the folder if it does not exists
        }
        xmf = fopen(extFilename.c_str(), "w");

        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"3.0\">\n");
        fprintf(xmf, "<!-- Kindly generated by poisson -->\n");
        fprintf(xmf, " <Domain>\n");
        fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");

        // the size of the grid has to be +1 since the 3DCoRectMesh is defined as vertex-centered
        fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n", topo->nglob(ax2) + 1, topo->nglob(ax1) + 1, topo->nglob(ax0) + 1);

        fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
        fprintf(xmf, "           <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "              %10.8f %10.8f %10.8f\n", 0.0, 0.0, 0.0);
        fprintf(xmf, "           </DataItem>\n");
        fprintf(xmf, "           <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "              %10.8f %10.8f %10.8f\n", 1.0, 1.0, 1.0);
        fprintf(xmf, "           </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n");

        if (!topo->isComplex()) {

            if (topo->lda() == 1){
                fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", attribute.c_str());
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0));
                fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), attribute.c_str());  //<-------------------------------
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
            } else {
                fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Cell\">\n", attribute.c_str());
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0), topo->lda());
                fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), attribute.c_str());  //<-------------------------------
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
            }
            
        } else {
            if (topo->lda() == 1){
                // real part
                string realname = attribute + "_real";
                fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", realname.c_str());
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0));
                fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), realname.c_str());  //<-------------------------------
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
                // imag part
                string imagname = attribute + "_imag";
                fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", imagname.c_str());
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0));
                fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), imagname.c_str());  //<-------------------------------
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
            } else {
                // real part
                string realname = attribute + "_real";
                fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Cell\">\n", realname.c_str());
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0), topo->lda());
                fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), realname.c_str());  //<-------------------------------
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
                // imag part
                string imagname = attribute + "_imag";
                fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Cell\">\n", imagname.c_str());
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0), topo->lda());
                fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), imagname.c_str());  //<-------------------------------
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
            }
        }

        //For Node centered:
        // fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", scalarName);
        // fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
        // fprintf(xmf, "        %s.h5:/data\n", filename); //<-------------------------------
        // fprintf(xmf, "       </DataItem>");
        // fprintf(xmf, "     </Attribute>");

        //For Vector:
        // fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Node\">\n", scalarName);
        // fprintf(xmf, "       <DataItem Dimensions=\"%d %d  3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
        // fprintf(xmf, "        xdmf2d.h5:/data"); //<-------------------------------
        // fprintf(xmf, "       </DataItem>");
        // fprintf(xmf, "     </Attribute>");
        //Also possible: Tensor (3x3), Tensor6 or even Matrix of NxM

        fprintf(xmf, "   </Grid>");
        fprintf(xmf, " </Domain>");
        fprintf(xmf, "</Xdmf>");
        fclose(xmf);
    }
    END_FUNC;
}

/**
 * @brief test the hdf5 dumping for real and complex numbers
 * 
 */
void hdf5_dumptest() {
    BEGIN_FUNC;

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int nglob[3] = {8, 8, 8};
    const int nproc[3] = {comm_size, 1, 1};

    //===========================================================================
    // real numbers
    Topology *topo    = new Topology(0, 1, nglob, nproc, false, NULL, 1, MPI_COMM_WORLD);
    const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};

    // we only allocate the real size = local size
    double *data = (double *)flups_malloc(sizeof(double *) * topo->locsize());

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                const size_t id = localIndex(0,i0, i1, i2, 0, 0, nmem,topo->nf());
                data[id + 0] = id;
            }
        }
    }
    // try the dump
    hdf5_dump(topo, "test_real", data);

    flups_free(data);
    delete (topo);

    //===========================================================================
    // create a real topology
    topo = new Topology(0, 1, nglob, nproc, true,NULL,1,MPI_COMM_WORLD);

    data = (double *)flups_malloc(sizeof(double *) * topo->locsize());

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                const size_t id = localIndex(0,i0, i1, i2, 0, 0, nmem,topo->nf());
                data[id + 0] = (double)id;
                data[id + 1] = 0.0;  //-data[id + 0];
            }
        }
    }
    // try the dump
    hdf5_dump(topo, "test_complex", data);
    flups_free(data);
    delete (topo);
    END_FUNC;
}