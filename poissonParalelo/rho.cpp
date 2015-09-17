 void H_rhoKernel (double *ne, double *ni, double *rho_h) {

    int size1 = J_X*J_Y*sizeof(cufftComplexDouble);

    cufftComplexDouble *rho_d;
    cudaMalloc(&rho_d, size1);
    cudaMemcpy(rho_d, rho_h, size1, cudaMemcpyHostToDevice);
    
    blockSize = 32; 
    dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE,1);
    dim3 dimGrid(ceil(float(J_Y/ blockSize),ceil(float( J_Y/ blockSize)), 1);
    D_rhoKernel<<< dimGrid, dimBlock >>>(ne, ni, rho_d)
    cudaDeviceSynchronize();

    cudaMemcpy(rho_h,rho_d, size1, cudaMemcpyDeviceToHost);
    cudaFree(rho_d);
  }

  __global__ 
  void D_rhoKernel(double *ne, double *ni, cufftComplexDouble *rho_d){

  int j = blockIdx.x*blockDim.x+threadIdx.x;
  int i = blockIdx.y*blockDim.y+threadIdx.y;
  int index;
  index = (i*J_Y+j);

  if((i<J_X)&&(j<J_Y)){
    rho[index].x= cte_rho* Factor_carga_e*(ni[index]- ne[index])/n_0;
    rho[index].y= 0.0;

  }

}
