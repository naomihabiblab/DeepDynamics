import torch
import torch.nn as nn


## Modeling
class ProbModelThin(nn.Module):

  def __init__(self, input_size=58):
    super(ProbModelThin, self).__init__()
    self.fc1 = nn.Linear(input_size, 99) # 58 bulk cell states
    self.fc2 = nn.Linear(99, 3)  # west prob, south prob, psuedotime
    self.relu = torch.relu
    self.softmax = torch.softmax

  def forward(self, x):
    x = self.relu(self.fc1(x))
    out = self.fc2(x)
    return out

  def predict(self, x):
    pred = self.forward(x)
    pred[:, :2] = self.softmax(pred[:, :2], dim=1)
    return pred
  

class ProbModel(nn.Module):

  def __init__(self, input_size=58):
    super(ProbModel, self).__init__()
    self.fc1 = nn.Linear(input_size, 90) # 58 cell states
    self.fc2 = nn.Linear(90, 140)
    self.fc3 = nn.Linear(140, 90)
    self.fc4 = nn.Linear(90, 3) # west prob, south prob, psuedotime
    self.relu = torch.relu
    self.softmax = torch.softmax

  def forward(self, x):
    x = self.relu(self.fc1(x))
    x = self.relu(self.fc2(x))
    x = self.relu(self.fc3(x))
    out = self.fc4(x)
    return out

  def predict(self, x):
    pred = self.forward(x)
    pred[:, :2] = self.softmax(pred[:, :2], dim=1)
    return pred
  

