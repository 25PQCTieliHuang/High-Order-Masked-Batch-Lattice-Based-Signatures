# High-Order Masked Batch Lattice-Based Signatures
This repository contains the source code for the paper **"High-Order Masked Batch Lattice-Based Signatures with Fast Masking Matrix-Vector Multiplication Gadget"**.

## Project Overview
We propose an efficient high-order masked batch signature scheme for lattice-based post-quantum cryptography, addressing the critical performance bottleneck of side-channel resistant implementations on resource-constrained embedded devices.

### Key Contributions
- **Fast Masking MVM Gadget**: A novel matrix-vector multiplication gadget based on the FastMMM algorithm, reducing polynomial multiplications by up to 33% compared to the state-of-the-art while maintaining provable t-NI security
- **Batch Signature Optimization**: An efficient batch signing method for the single-user multi-message setting, with offline public-key preprocessing and negligible memory overhead
- **Generalizable Design**: Ported to both Dilithium (NIST PQC standard) and Raccoon signature schemes
- **Comprehensive Evaluation**: Full C implementation with detailed performance benchmarks across different masking orders (1st to 31st order)

### Performance Highlights
- **Dilithium**: Up to 21.26% speedup for the masking MVM operation
- **Raccoon**: Up to 35.16% performance gain for batch signing with 20 messages
- **Memory Efficiency**: Memory overhead remains nearly identical to single-message signing

## Compilation Environment
### System Requirements
- **Operating System**: Linux (Ubuntu 22.04 LTS or later recommended)
- **Processor**: x86_64 architecture (tested on Intel Core i7-12800HX)
- **Memory**: Minimum 4GB RAM (8GB recommended for high-order masking)

### Software Dependencies
| Dependency | Version |
|------------|---------|
| GCC | 11.4.0 or higher |
| CMake | 3.16 or higher |
| Make | 4.3 or higher |


## License
This project is licensed under the MIT License - see the LICENSE file for details.
