# Predict Storage Densities For All 231 Billion
This is a Rust program that uses the linear regression models (storage-24_model.csv and storage-57_model.csv) to estimate the storage density for all 230 billion molecules. 
The program writes genes for molecules with a storage density > 0.30 kJ/g to the file all_230B_data.dat.

But first, compile the program:
```
cargo build --release
```
Then run the new executable and wait (~12 hours).