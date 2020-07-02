task qc{
	input {	
		File    infile
		String  filetype
		String	project_name
		String?  cpu='2'
		String?  memory='4G'
		Int     min_gene
		Int     min_cell
		Int     genes 
		Int     per_mt
	}
	command <<<
		set -e
		set -o pipefail
		python << code
		import scanpy as sc
		filetype = "~{default="csv" filetype}"
		infile   = "~{infile}"
		if filetype=='csv':
			adata = sc.read_csv(infile, delimiter=",", first_column_names=True).T
		elif filetype=='tsv':
			adata = sc.read_csv(infile, delimiter="\t", first_column_names=True).T
		elif filetype=='10x_mtx':
			adata = sc.read_csv(infile)
		sc.settings.autosave = True
		#filter and QC
		sc.pl.highest_expr_genes(adata, n_top=20, save=".png")
		sc.pp.filter_cells(adata, min_genes=~{min_gene})
		sc.pp.filter_genes(adata, min_cells=~{min_cell})
		adata.var['mt'] = adata.var_names.str.startswith('MT-')
		sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
		sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save=".png")
		sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save=".png")
		sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save=".png")
		adata = adata[adata.obs.n_genes < ~{genes}, :]
		adata = adata[adata.obs.pct_counts_mt < ~{per_mt}, :]
		adata.write("~{project_name}.h5ad", compression="gzip")	
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
    	File outputfile = "${project_name}.h5ad"
 	Array[File] pngfile    = glob("figures/*png")
  }
}

