
# Merge_PDF_first_pages in Batch_LCModel_Toolkit
# Author: Kelley Swanberg (Columbia University, 2023)
# Written for R Version 4.0.5
# Merge the first page of every PDF in a chosen directory into a single master PDF. 
# Enables easy batch perusal of spectral fit visuals in LCModel output .ps files converted to PDF 

library('staplr')
library('stringr')

# Allow user to select directory of PDFs 
working_directory = choose.dir(default = "", caption = "Select folder")
setwd(working_directory); 

pdf_files = list.files(path = working_directory, pattern = "*.pdf"); 
pdf_files_to_merge = pdf_files; 
length_pdf_files = length(pdf_files); 

for (ii in 1:length_pdf_files){
  file_path = pdf_files[ii]
  saving_path = str_replace(file_path, '.pdf', '_pg1.pdf')
  pdf_files_to_merge[ii] = saving_path; 
  staplr::select_pages(selpages = 1, input_filepath = file_path, output_filepath = saving_path)
}

# Merge all extracted first pages 
staplr::staple_pdf(input_files=pdf_files_to_merge, output_filepath='merged_pdfs.pdf')