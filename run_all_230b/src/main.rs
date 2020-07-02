use std::fs;
use std::io::prelude::*;
use std::io::LineWriter;

extern crate ndarray;
use ndarray::{Array, Array1};

fn read_coefficients(filename: String) -> (f64, Array1<f64>) {
    // Read coefficient of linear regression model.
    // First value are the intersept, the remaning corresponds to the regression coefs.

    let coefs = fs::read_to_string(filename)
        .expect("Error reading file");

    let mut collect_coefs: Vec<f64> = Vec::with_capacity(287);
    let mut intersept: f64 = 0.;
    for (i, coef) in coefs.split_whitespace().enumerate() {
        if i == 0 {
            intersept += coef.parse::<f64>().unwrap();
            continue
        }
        collect_coefs.push(coef.parse::<f64>().unwrap());
    }
    (intersept, Array::from_vec(collect_coefs))
}

fn make_gene(idxs: &[usize; 7]) -> Array1<f64> {
    // Make gene - return 1x287 matrix - which corresponds to the flattened 7x41.

    let mut gene = Array1::<f64>::zeros(287);
    for (i, sub_idx)  in idxs.iter().enumerate() {
        if sub_idx == &0 {
            continue;
        }
        gene[[i*41 + sub_idx-1]] = 1.;
    }
    gene 
}

fn compute_storage_energy(idx: &[usize; 7], coefs: &Array1<f64>, intersept: &f64) -> f64 {
    // Compute storage density of gene.
    let gene = make_gene(idx);
    gene.dot(coefs) + intersept
}

fn write_genes(genes: &Vec<[usize; 7]>, densities: &Vec<f64>, outfile: String) {
    // Write genes to file. 

    let file = fs::File::create(outfile).expect("can't open write file");
    let mut file = LineWriter::new(file);

    for (gene, density) in genes.iter().zip(densities.iter()) {
        let mut gene_data: String = gene.into_iter().map(|i| format!("{}-",i).to_string())
                                                    .collect::<String>();
        gene_data.pop();
        gene_data.push(',');
    
        gene_data.push_str(&density.to_string());
        gene_data.push_str("\n");
        file.write_all(gene_data.as_bytes()).expect("can't write gene")
    }
}


fn main() {

    // Load all coefficients.
    let (intersept24, coefs24) = read_coefficients("storage-24_model.csv".to_string());
    let (intersept57, coefs57) = read_coefficients("storage-57_model.csv".to_string());
    
    assert!(intersept24 == 0.241495481106874, "intersept of 1-4 sub not OK!");
    assert!(intersept57 == 0.306111038383119, "intersept of 5-7 sub not OK!");

    // Check that some of the densities are OK.
    let parrent_gene: [usize; 7] = [41, 0, 0, 0, 0 ,0 ,0];
    assert!(compute_storage_energy(&parrent_gene, &coefs24, &intersept24) == 0.259466163684428, "Parent energy not OK");
    let test_gene1: [usize; 7] = [13, 0, 0, 40, 0 ,0 ,1];
    assert!(compute_storage_energy(&test_gene1, &coefs24, &intersept24) == 0.212504288135473, "test1 energy not OK");
    let test_gene2: [usize; 7] = [13, 0, 21, 40, 0 , 6, 1];
    assert!(compute_storage_energy(&test_gene2, &coefs57, &intersept57) == 0.219730651174742, "test2 energy not OK");
    
    // Run all 230B and collect genes, and energies
    let mut collect_gene: Vec<[usize; 7]> = Vec::with_capacity(100_000_000);
    let mut collect_density: Vec<f64> = Vec::with_capacity(100_000_000);

    let mut counter: i64 = 0;
    for i in 0..42 {
        for j in 0..42 {
            for k in 0..42 { 
                for l in 0..42 {
                    for m in 0..42 {
                        for n in 0..42 {
                            for o in 0..42 {
                                
                                let gene: [usize; 7] = [i, j, k, l, m ,n , o];
                                let mut density: f64 = 0.0;

                                if 7 - &gene.iter().filter(|&n| *n == 0).count() <= 4 {
                                    density += compute_storage_energy(&gene, &coefs24, &intersept24);
                                } else {
                                    density += compute_storage_energy(&gene, &coefs57, &intersept57);
                                }

                                if density >= 0.30 {
                                    collect_gene.push(gene);
                                    collect_density.push(density);
                                }

                                counter += 1;
                                if counter % 1_000_000_000 == 0 {
                                    println!("iter {} of {}.", counter, 42f64.powf(7.))
                                }

                            }
                        }
                    }
                }
            }
        }
    }

    // Write genes with storage above 0.25 kJ/g.
    write_genes(&collect_gene, &collect_density, "all_230B_data.dat".to_string())
}
