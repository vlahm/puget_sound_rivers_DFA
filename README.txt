When using any of my R files, you'll want to collapse folds. This makes them much easier to navigate. use Alt+O on Windows/Linux or Cmd+Opt+O on Mac
all subfolders contain their own READMEs, but this is a good summary.

00_tmb_uncor_Rmat.R contains Tim Cline's TMB script and associated functions for Dynamic Factor Analysis. 
	It is called by 02_ and 03_, and need not be touched.
01_preprocessing is for structuring the chemical/physical stream data, and should probably be in the associated data folder, but here it is.
	No need to change it either.
02_testing_or_evaluation.R and 03_model_fitting.R are set up similarly, but the former is used for tweaking and testing models before they are
	used in the giant model_fitting loop contained within 03_. 03 will save copies of all model objects as .rds files inside model_objects.
	pdfs of trends, loadings, fits, and regressions will be in model_outputs.
the results of the best model-fitting round are in round_6_TeTuSu_UNSCALED.
After the best model has been selected, it can be further explored in 02_ (bottom sections). You'll want to structure the top of 02_ as if you're
	about to test that model (i.e. choose the same number of trends, same covariates, etc.), but then skip section 4.
saved_structures can be used to store secondary model objects so that they don't have to be regenerated during exploration.
	You can load objects from here in section 4.2
transform_testing was used to determine that backtransforming covariate effect sizes into their original scales/shapes after running the model is not
	possible. This pertains to both power transformation and box-cox transformation. Log too, I presume. Instead, effect sizes are simply reported
	as log(change_response)/change_covariate
NHDPlus_variables_and_data_sources contains extra info on watershed variables and their data sources
