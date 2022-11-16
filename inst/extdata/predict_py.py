
from tensorflow.keras.models import load_model
import pandas as pd
import numpy as np

def predict_py(path_score,Path):
  data2 = path_score
  model = load_model(Path+"/TCfinder.hdf5")
  predict = model.predict(data2)
  return predict
