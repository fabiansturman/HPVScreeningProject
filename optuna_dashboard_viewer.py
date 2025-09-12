from optuna_dashboard import run_server
import optuna

name = "C:\\Users\\fabia\\Documents\\Uni\\DPhil\\SophieHPV\\HPVsim_FABIAN\\CalibrationRawResults\\d2Cal_8Sep25_2"


run_server(f'sqlite:///{name}.db', host="127.0.0.1", port = 8084)