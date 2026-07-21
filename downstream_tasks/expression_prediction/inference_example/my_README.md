Kill all processes
pkill -u dpanc -f ipykernel
pkill -u dpanc -f torch/_inductor/compile_worker

Use this notebook for inference, there is batch size which you could adjust

/home/dpanc/benchmarking/GENA_LM/GENA_LM_expression_branch/downstream_tasks/expression_prediction/inference_example/inference_polina_batch.ipynb

Then lets try notebook of Egor
upload it form git
wget -O expression_inference_example_api.ipynb \
https://raw.githubusercontent.com/AIRI-Institute/GENA_LM/api/downstream_tasks/expression_prediction/api/notebooks/expression_inference_example.ipynb


