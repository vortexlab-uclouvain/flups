#include <mpi.h>

#include "gtest/gtest.h"

int main(int argc, char** argv) {
    /* Initialize the infrastructure necessary for communication */
    MPI_Init(&argc, &argv);
    //-------------------------------------------------------------------------
    int err;
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        for (int ia = 0; ia < argc; ia++) {
            int comp;
            comp = strcmp(argv[ia], "--gtest_output=xml:report_test.xml");
            if (comp == 0) {
                sprintf(argv[ia], "--gtest_output=xml:report_test_rank%d.xml", rank);
            }

            comp = strcmp(argv[ia], "--gtest_output=xml:report_valid.xml");
            if (comp == 0) {
                sprintf(argv[ia], "--gtest_output=xml:report_valid_rank%d.xml", rank);
            }
        }
        printf("Starting the tests...  \n");
        // init google and remove Google Test arguments
        ::testing::InitGoogleTest(&argc, argv);
        err = RUN_ALL_TESTS();
    }
    //-------------------------------------------------------------------------
    /* Tear down the communication infrastructure */
    MPI_Finalize();
    return err;
}