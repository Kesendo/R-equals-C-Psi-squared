@echo off
rem Detached driver for the N=9 layer discharge: sector E, then O (checkpointed, resumable).
cd /d "D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
python -u simulations\_n9_layer_discharge.py E > simulations\_n9_discharge_E_log.txt 2>&1
python -u simulations\_n9_layer_discharge.py O > simulations\_n9_discharge_O_log.txt 2>&1
echo done > simulations\_n9_discharge_DONE.txt
