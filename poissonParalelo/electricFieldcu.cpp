 void H_electricFiled (double *phi_h, double *ex_h double *ey_h, double hx) {

    int size1 = J_X*J_Y*sizeof(double);

    double *ex_d, *ey_d, *phi_d;
    cudaMalloc(&ex_d, size1);
    cudaMalloc(&ey_d, size1);
    cudaMalloc(&phi_d, size1);
    cudaMemcpy(phi_d, phi_h, size1, cudaMemcpyHostToDevice);
    cudaMemset(ex_d, 0, size1);
    cudaMemset(ey_d, 0, size1);

    dim3 dimBlock(BLOCK_SIZE,1,1);
    dim3 dimGrid3(ceil(float(J_X)/ blockSize), 1, 1);
    D_electricFiled<<< dimGrid3, dimBlock >>> (phi_d, ex_d, ey_d, hx)
    cudaDeviceSynchronize();

    cudaMemcpy(ex_h,ex_d, size1, cudaMemcpyDeviceToHost);
    cudaMemcpy(ey_h,ey_d, size1, cudaMemcpyDeviceToHost);
    cudaFree(ex_d);
    cudaFree(ey_d);
    cudaFree(phi_d);
  }

  __global__
  void D_electricFiled (double *phi_d, double *ex_d, double *ey_d, double hx) {

    int j = (blockIdx.x * blockDim.x + threadIdx.x)+1;

     if((j<J_X-1)){

      for (int k=1;k<J_Y-1;k++)
      {
        ex_d[j*J_Y+k]=(phi_d[(j-1)*J_Y+k]-phi_d[(j+1)*J_Y+k])/(2.*hx);
        ey_d[j*J_Y+k]=(phi_d[j*J_Y+(k-1)]-phi_d[j*J_Y+(k+1)])/(2.*hx);

        ex_d[0*J_Y+k] = 0.0;  //Cero en las fronteras X
        ey_d[0*J_Y+k] = 0.0;
        ex_d[(J_X-1)*J_Y+k] = 0.0;
        ey_d[(J_X-1)*J_Y+k] = 0.0;
      }

      ex_d[j*J_Y+0]=(phi_d[(j-1)*J_Y+0]-phi_d[((j+1)*J_Y+0)])/(2.*hx);
      ey_d[j*J_Y+0]=(phi_d[j*J_Y+(J_Y-1)]-phi_d[j*J_Y+1])/(2.*hx);


      ex_d[j*J_Y+(J_Y-1)]=(phi_d[(j-1)*J_Y+(J_Y-1)]-phi_d[(j+1)*J_Y+(J_Y-1)])/(2.*hx);
      ey_d[j*J_Y+(J_Y-1)]=(phi_d[j*J_Y+(J_Y-2)]-phi_d[j*J_Y+0])/(2.*hx);
  }
  }