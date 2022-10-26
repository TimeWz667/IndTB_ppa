from ppa.dy.base import AbsModel

__author__ = 'Chu-Chang Ku'
__all__ = ['I']


class I:
    Asym = 0
    Sym = 1
    ExCsPub = 2
    ExCsEng = 3
    ExCsPri = 4
    TxPub = 5
    TxEng = 6
    TxPri = 7
    Cured = 8
    SelfCured = 9
    LTFU = 10
    DeadU = 11
    DeadT = 12

    N_States = 13

