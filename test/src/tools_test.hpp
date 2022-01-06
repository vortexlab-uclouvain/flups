#ifndef _SRC_TOOLS_TEST_
#define _SRC_TOOLS_TEST_

static const FLUPS_BoundaryType FLUPS_bcs[3] = {EVEN, ODD, UNB};

static void GetDirIdxFromId(int glob_id, int bcidx[3]){
    const int ncomb = 10; // number of different combination of boundary condition for 1 dimension
    const int n0    = ncomb;
    const int n01   = ncomb * n0;
    const int n012  = ncomb * n01;

    bcidx[2] = glob_id / n01;
    bcidx[1] = (glob_id % n01) / n0;
    bcidx[0] = (glob_id % n01) % n0;
}

static void GetBcIdxFromId(int glob_id, FLUPS_BoundaryType bcs[2]){
    if(glob_id == 9){
        bcs[0] = PER;
        bcs[1] = PER;
    }else{
        int i0 = glob_id%3; 
        int i1 = glob_id/3;
        bcs[0] = FLUPS_bcs[i0]; 
        bcs[1] = FLUPS_bcs[i1];
    }
}

static void InitBcFromId(int glob_id, FLUPS_BoundaryType* bc[3][2]){
    int dir_idx[3];
    GetDirIdxFromId(glob_id, dir_idx);    

    // For the moment, we test only scalar field
    for(int dir = 0; dir < 3; dir ++){
        FLUPS_BoundaryType bcs[2];
        GetBcIdxFromId(dir_idx[dir], bcs);
        bc[dir][0][0] = bcs[0]; //FLUPS_bcs[bcidx[0]];
        bc[dir][1][0] = bcs[1]; //FLUPS_bcs[bcidx[1]];    
        
    }
}


#endif //_SRC_TOOLS_TEST_

