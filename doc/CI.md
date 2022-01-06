# Automatic testing using Gitlab Continuous Integration Tool

Most of the information have been found in the Gitlab documentation ([here](https://docs.gitlab.com/ee/ci/)).  

The Continuous Integration of this project relies on 2 tools: the Gitlab CI and the [Google Test framework](https://docs.gitlab.com/ee/ci/). The latter organizes individual tests into a series of test suites and produces test report which can be processed by gitlab. 

## Addition of test in Gitlab 

The configuration options for Gitlab can be found in the `.gitlab-ci.yml` [file](../.gitlab-ci.yml). 
Gitlab runs the tests using a Gitlab runner and a docker image. The runner build and launch the docker image, copy the source of the current commit of the project inside the workspace and perform the tests. Therefore, when creating new test, you first need to __make sure that the docker image you are using contains all the libraries needed for the tests.__

In this project, on jan. 6 2022, the image is `immc/flups-valid:0.5` and it has the following libraries: 
- fftw
- hdf5
- metis
- google test 
The docker file used to generate the image can be found [here](../docker/DockerFile)

Then, you can start the configuration of your CI. The jobs (fundamental element of the `.yml` file, giving the instruction of for a specific task) can be grouped in different stage, which describe the sequential execution of jobs. In this project, we have choosen 4 different stages (`build`, `test`, `vector` and `validation`), each of them corresponding to a certain type of tests. 


## Use the cluster to perform some test 

Sometimes, it may be usefull to run some tests on cluster. There are some tricks to do that. 

First, you will need a username and a password to access the selected cluster. Once you have it, you need to store them on Gitlab (Settings->CI->Variables). For security reasons, those variables should be stored in a protected state. Thoses variable can be used directly in the `.gitlab-ci.yml`. 

Then, you need to add the fingerprint of the cluster in the knownhost of your image, and the key in you `.ssh/config`. This will enable you to access the cluster from your image. Then, you 'll only need to implement the jobs you want to perform on the cluster. 