import torch
import torch.nn as nn
import torch.nn.functional as F


class ProbLoss(nn.Module):
    def __init__(self, alpha=1, beta=1, gamma=1):
        super(ProbLoss, self).__init__()
        self.__alpha = alpha
        self.__beta = beta
        self.__gamma = gamma
        self.__kl = nn.KLDivLoss(reduction="batchmean")
        self.__ce = nn.CrossEntropyLoss(reduction='mean')
        self.__time_loss = nn.MSELoss()

    def forward(self, output, target):
        log_output = torch.log_softmax(output, dim=1)
        kl_loss = self.__kl(log_output[:, :2], target[:, :2])
        ce_loss = self.__ce(output[:, :2], target[:, :2])
        time_loss = self.__time_loss(output[:, 2], target[:, 2])
        loss = self.__alpha * ce_loss + self.__beta * kl_loss + self.__gamma * time_loss
        return loss, kl_loss, ce_loss, time_loss

