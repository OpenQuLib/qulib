  Library of Qulib

  What is it?
  -----------
     o  The library of Qulib is a mpi parallelization function library
        for simulation  quantum computing, which is mainly used for 
        multi-machine parallelization classical simulation of quantum 
        computing. By adopting a distributed storage method, the  
        problem that the occupied memory space in the classical computer 
        simulation quantum computing process increases exponentially 
        with the increase of the number of analog quantum bits is 
        optimized to some extent. The main implementation of the mpi 
        parallelization is as follows:
   
     o  Initialization: Initialize the quantum state of the fixed number 
        of qubits and store the probability amplitude distribution on 
        each node.

     o  Print: Outputs the probability amplitude of the quantum state 
        stored in the quantum   register.

     o  Delete: Empty the quantum state in the quantum register and 
        reclaim the memory.

     o  MPI parallelization calculation of correlation function: 
        quantum_qft function, quantum_gate1 function, quantum_sigma_x, 
        quantum_sigma_y, quantum_sigma_z operation, quantum_toffoli, 
        quantum_cnot gate-operation. Quantum_bmeasure,quantum_qft_state_collapse 
        function.
 

  Contacts
  --------

     o If you want to be informed about new code releases, bug fixes,
       security fixes, general news and information about the library
       subscribe to the qulib-announce mailing address as described under
       <openqulib@gmail.com>

     o If you want freely available support for running the library,please see the
       resources at <http://qulib.org>

     o If you have a concrete bug report for Qulib please see the instructions
       for bug reporting at <http://qulib.org/bug_report.html>


