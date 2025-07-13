
import pandas as pd
import torch
from torch.utils.data import Dataset
from sklearn.model_selection import train_test_split

PATH_TO_DATA = '/ems/elsc-labs/habib-n/yuval.rom/500/DeepDynamics/prediction/data/shared_bulk_data_mask.csv'
PATH_TO_TARGET = '/ems/elsc-labs/habib-n/yuval.rom/500/DeepDynamics/prediction/data/y.csv'

TRAIN_PATH = 'data/train_data.pt'
TEST_PATH = 'data/test_data.pt'


class SCData(Dataset):
  def __init__(self, X, y):
    self.X = X
    self.y = y

  def __getitem__(self, index):
    return self.X[index], self.y[index]

  def __len__(self):
    return len(self.X)


if __name__ == '__main__':
  # Load the data
    shared_bulk_data_mask = pd.read_csv(PATH_TO_DATA, index_col=0)
    y_mask = pd.read_csv(PATH_TO_TARGET, index_col=0)
    input_size = shared_bulk_data_mask.shape[1]
    shared_bulk_data_mask.shape
    X_train, X_test, y_train, y_test = train_test_split(shared_bulk_data_mask, y_mask, test_size=0.25, random_state=42)

    X_train_tensor = torch.tensor(X_train.values, dtype=torch.float32)
    y_train_tensor = torch.tensor(y_train.values, dtype=torch.float32)
    X_test_tensor = torch.tensor(X_test.values, dtype=torch.float32)
    y_test_tensor = torch.tensor(y_test.values, dtype=torch.float32)

    train_data = SCData(X_train_tensor, y_train_tensor)
    test_data = SCData(X_test_tensor, y_test_tensor)

    # Save the datasets
    torch.save(train_data, TRAIN_PATH)
    torch.save(test_data, TEST_PATH)