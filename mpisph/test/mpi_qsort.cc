TEST(){
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  srand(time(NULL)*rank);

  // Generating the particles randomly on each process
  size_t nparticles = 2000;
  size_t nparticlesperproc = nparticles/size; 
  double maxbound = 1.0; // Particles positions between [0,1]
  // Adjust for last one 
  if(rank == size-1){
    nparticlesperproc = (nparticles-nparticlesperproc*(size-1));
  }
  if(rank==0){
    std::cout<<"Generating "<<nparticles<<std::endl;
  }
  std::cout<<"Rank "<<rank<<nparticlesperproc<<" particles"<<std::endl;

  // Range to compute the keys 
  
	std::vector<std::pair<entity_key_t,body>> bodies(nparticlesperproc);
  // Create the bodies and keys 
  for(size_t i=0;i<nparticlesperproc;++i){
    // Random x, y and z
    bodies[i].second.setPosition(
      (double)rand()/(double)RAND_MAX*(maxbound),
      (double)rand()/(double)RAND_MAX*(maxbound),
      (double)rand()/(double)RAND_MAX*(maxbound)
    );
 
    // Compute the key 
    bodies[i].first = entity_key_t(range,bodies[i].second.getPosition());
  }

  // Gather all the particles everywhere and sort locally 
  std::vector<std::pair<entity_key_t,body>> checking(nparticles);
  MPI_Allgather(
    &bodies[0],
    nparticlesperproc*sizeof(std::pair<entity_key_t,body>),
    MPI_BYTE,
    &checking[0],
    nparticlesperproc*sizeof(std::pair<entity_key_t,body>),
    MPI_BYTE,
    MPI_COMM_WORLD);

  // Sort it locally base on the keys 
  std::sort(checking.begin(),checking.end(),
      [](auto& left, auto& right){
        return left.first < right.first;
      });

  // Extract the subset of this process 
  std::vector<std::pair<entity_key_t,body>> my_checking(
    checking.begin()+rank*(nparticles/size),
    checking.begin()+rank*nparticlesperproc+nparticlesperproc);

	// Use the mpi_qsort
  mpi_qsort(bodies,nparticles); 
	
	// Compare the results with all processes particles subset 
  CINCH_ASSERT(my_checking == bodies);
}