import numpy as np
import torch
import torch.optim as optim
from torch.utils.data import DataLoader

from prediction.loss import *
from prediction.model import *
from prediction.data_structures import *


TRAIN_PATH = 'data/train_data.pt'
TEST_PATH = 'data/test_data.pt'

OUTPUT_PATH = 'outputs/'
WEIGHTS_PATH = 'weights/'


def eval_model(model, optimizer, test_loader, loss_func):
    total_loss = 0
    kl_total_loss, ce_total_loss, time_total_loss = 0, 0, 0
    N = len(test_loader.dataset) / test_loader.batch_size
    model.eval()
    with torch.no_grad():
      for batch_X, batch_y in test_loader:
        output = model(batch_X)
        loss, kl_loss, ce_loss, time_loss = loss_func(output, batch_y)

        total_loss += loss
        kl_total_loss += kl_loss
        ce_total_loss += ce_loss
        time_total_loss += time_loss

    model.train()
    return total_loss.item() / N, kl_total_loss.item() / N, ce_total_loss.item() / N, time_loss.item() / N


def train_one_epoch(model, optimizer, train_loader, loss_func, config):
    epoch_loss, epoch_kl_loss, epoch_ce_loss, epoch_time_loss = 0, 0, 0, 0
    N = len(train_loader.dataset) / train_loader.batch_size
    for batch_X, batch_y in train_loader:
      output = model(batch_X)
      loss, kl_loss, ce_loss, time_loss = loss_func(output, batch_y)

      if config:
        l1_reg = sum(param.abs().sum() for param in model.parameters())
        l2_reg = sum(param.pow(2).sum() for param in model.parameters())
        loss += config['l1_lambda'] * l1_reg + config['l2_lambda'] * l2_reg

      optimizer.zero_grad()
      loss.backward()
      optimizer.step()

      epoch_loss += loss
      epoch_kl_loss += kl_loss
      epoch_ce_loss += ce_loss
      epoch_time_loss += time_loss

    return epoch_loss.item() / N, epoch_kl_loss.item() / N, epoch_ce_loss.item() / N, epoch_time_loss.item() / N


def train_model(model, optimizer, train_loader, test_loader, epochs_num,
                loss_func, print_progress=True, early_stopping=True, patience=10, config=None):
  train_loss, test_loss = [], []
  kl_train_loss_list, ce_train_loss_list, time_train_loss_list = [], [], []
  kl_test_loss_list, ce_test_loss_list, time_test_loss_list = [], [], []
  best_loss, counter = np.inf, 0

  for epoch in range(epochs_num):
    epoch_loss, epoch_kl_loss, epoch_ce_loss, time_train_loss = train_one_epoch(model,
                                                               optimizer,
                                                               train_loader,
                                                               loss_func,
                                                              config)

    eval_loss, kl_eval_loss, ce_eval_loss, time_eval_loss = eval_model(model, optimizer,
                                                       test_loader, loss_func)
    train_loss.append(epoch_loss)
    test_loss.append(eval_loss)
    kl_train_loss_list.append(epoch_kl_loss)
    kl_test_loss_list.append(kl_eval_loss)
    ce_train_loss_list.append(epoch_ce_loss)
    ce_test_loss_list.append(ce_eval_loss)
    time_train_loss_list.append(time_train_loss)
    time_test_loss_list.append(time_eval_loss)
    # early stopping
    if eval_loss < best_loss:
        best_loss = eval_loss
        counter = 0
    else:
        counter += 1
    if counter >= patience and early_stopping:
        print(f'Early stopping at epoch {epoch + 1}')
        print(f"Epoch {epoch + 1}: Train Loss = {epoch_loss:.3f}, Test Loss = {eval_loss:.3f}")
        break
    # print progress
    if (epoch + 1) % 10 == 0 and print_progress:
      print(f"Epoch {epoch + 1}: Train Loss = {epoch_loss:.3f}, Test Loss = {eval_loss:.3f}")

  return (train_loss, test_loss, kl_train_loss_list, ce_train_loss_list,
          kl_test_loss_list, ce_test_loss_list,
          time_train_loss_list, time_test_loss_list)


if __name__ == "__main__":
   
    # Load the data
    train_data = torch.load(TRAIN_PATH, weights_only=False)
    test_data = torch.load(TEST_PATH, weights_only=False)

    #### Model settings
    batch_size = 10
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=batch_size, shuffle=False)

    config = {"l1_lambda": 5.541763126909059e-05,
            "l2_lambda": 0,
            }
    # config = {"l1_lambda": 0,
    #         "l2_lambda": 0.,
    #         }
    
    model = ProbModel()
    loss_func = ProbLoss()
    optimizer = optim.Adam(model.parameters(), lr=1e-5)
    epochs = 2000
    train_loss, test_loss, kl_train_loss_list, ce_train_loss_list, kl_test_loss_list, ce_test_loss_list, time_train_loss_list, time_test_loss_list = train_model(model, optimizer, train_loader, test_loader, epochs, loss_func, patience=10, config=config)
    
    # save model and results
    torch.save(model.state_dict(), WEIGHTS_PATH + 'model.pth')
    np.save(OUTPUT_PATH + 'train_loss.npy', np.array(train_loss))
    np.save(OUTPUT_PATH + 'test_loss.npy', np.array(test_loss))
    np.save(OUTPUT_PATH + 'kl_train_loss.npy', np.array(kl_train_loss_list))
    np.save(OUTPUT_PATH + 'ce_train_loss.npy', np.array(ce_train_loss_list))
    np.save(OUTPUT_PATH + 'kl_test_loss.npy', np.array(kl_test_loss_list))
    np.save(OUTPUT_PATH + 'ce_test_loss.npy', np.array(ce_test_loss_list))
    np.save(OUTPUT_PATH + 'time_train_loss.npy', np.array(time_train_loss_list))
    np.save(OUTPUT_PATH + 'time_test_loss.npy', np.array(time_test_loss_list))


    # model_thin = ProbModelThin()
    # optimizer = optim.Adam(model_thin.parameters(), lr=1e-4)
    # epochs = 2000
    # train_loss_thin, test_loss_thin, kl_train_loss_list_thin, ce_train_loss_list_thin, kl_test_loss_list_thin, ce_test_loss_list_thin, time_train_loss_list_thin, time_test_loss_list_thin = train_model(model_thin, optimizer, train_loader, test_loader, epochs, loss_func, patience=10, config=config)

    # torch.save(model_thin.state_dict(), WEIGHTS_PATH + 'model_thin.pth')
    # np.save(OUTPUT_PATH + 'train_loss_thin.npy', np.array(train_loss_thin))
    # np.save(OUTPUT_PATH + 'test_loss_thin.npy', np.array(test_loss_thin))
    # np.save(OUTPUT_PATH + 'kl_train_loss_thin.npy', np.array(kl_train_loss_list_thin))
    # np.save(OUTPUT_PATH + 'ce_train_loss_thin.npy', np.array(ce_train_loss_list_thin))
    # np.save(OUTPUT_PATH + 'kl_test_loss_thin.npy', np.array(kl_test_loss_list_thin))
    # np.save(OUTPUT_PATH + 'ce_test_loss_thin.npy', np.array(ce_test_loss_list_thin))
    # np.save(OUTPUT_PATH + 'time_train_loss_thin.npy', np.array(time_train_loss_list_thin))
    # np.save(OUTPUT_PATH + 'time_test_loss_thin.npy', np.array(time_test_loss_list_thin))
