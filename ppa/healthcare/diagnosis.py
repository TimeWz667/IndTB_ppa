from collections import Counter
import numpy.random as rd


__author__ = 'Chu-Chang Ku'
__all__ = ['Algorithm', 'Test', 'Specimen', 'Results']


class Results:
    def __init__(self):
        self.TruePos = 0
        self.FalseNeg = 0
        self.FalseNegPreDx = 0
        self.FalsePos = 0
        self.TrueNeg = 0
        self.TrueNegPreDx = 0
        self.Counts = Counter()

    def __getitem__(self, key):
        return self.Counts[key]

    def count(self, k, n):
        self.Counts[k] += n

    def push(self, tp=0, fn=0, fn_pd=0, fp=0, tn=0, tn_pd=0):
        self.TruePos += tp
        self.FalseNeg += fn
        self.FalseNegPreDx += fn_pd
        self.FalsePos += fp
        self.TrueNeg += tn
        self.TrueNegPreDx += tn_pd

    @property
    def ppv(self):
        try:
            return self.TruePos / (self.FalsePos + self.TruePos)
        except ZeroDivisionError:
            return 0

    @property
    def neg_predx(self):
        return (self.FalseNegPreDx + self.TrueNegPreDx) / (self.FalseNeg + self.TrueNeg)

    def print(self, detail=False):
        print(f'-- TP/FN (FN pre-dx) = {self.TruePos:.0f}/{self.FalseNeg:.0f} ({self.FalseNegPreDx:.0f})')
        print(f'-- FP/TN (TN pre-dx)= {self.FalsePos:.0f}/{self.TrueNeg:.0f} ({self.TrueNegPreDx:.0f})')
        print(f'-- PPV = {self.ppv:.2%}')
        print(f'-- Negative due to pre-dx LTFU = {self.neg_predx:.2%}')

        if detail:
            for k, v in self.Counts.items():
                print(f'-- {k}: {v}')

    def __add__(self, other):
        res = Results()
        res.push(tp=self.TruePos, fn=self.FalseNeg, fn_pd=self.FalseNegPreDx,
                 fp=self.FalsePos, tn=self.TrueNeg, tn_pd=self.TrueNegPreDx)
        res.Counts.update(self.Counts)
        res.push(tp=other.TruePos, fn=other.FalseNeg, fn_pd=other.FalseNegPreDx,
                 fp=other.FalsePos, tn=other.TrueNeg, tn_pd=other.TrueNegPreDx)
        res.Counts.update(other.Counts)
        return res


class Specimen:
    def __init__(self, name, loss):
        self.Name = name
        self.Loss = loss

    def collect(self, n):
        return (n * (1 - self.Loss), n * self.Loss), (f'N_collect_{self.Name}', n)

    def collect_sto(self, n):
        loss = rd.binomial(n, p=self.Loss)
        return (n - loss, loss), (f'N_collect_{self.Name}', n)


class Test:
    def __init__(self, name, sens, spec, specimen=None):
        self.Name = name
        self.Sens = sens
        self.Spec = spec
        self.Specimen = specimen

    def test_tb(self, n):
        return (n * self.Sens, n * (1 - self.Sens)), (f'N_test_{self.Name}', n)

    def test_nontb(self, n):
        return (n * (1 - self.Spec), n * self.Spec), (f'N_test_{self.Name}', n)

    def test_tb_sto(self, n):
        pos = rd.binomial(n, p=self.Sens)
        return (pos, n - pos), (f'N_test_{self.Name}', n)

    def test_nontb_sto(self, n):
        neg = rd.binomial(n, p=self.Spec)
        return (n - neg, neg), (f'N_test_{self.Name}', n)


class Algorithm:
    def __init__(self, key, ssm: Test = None, xpert: Test = None, cdx: Test = None):
        self.Key = key
        self.SSM = ssm
        self.Xpert = xpert
        self.CDx = cdx

    def dx(self, n_tb=0, n_nontb=0):
        res = Results()

        tp, fn, fnl = 0, n_tb, 0
        fp, tn, tnl = 0, n_nontb, 0

        specimen = None
        if self.SSM is not None:
            specimen = self.SSM.Specimen
            if specimen is not None:
                (fn, fnl1), stats = specimen.collect(fn)
                fnl += fnl1
                res.count(stats[0], stats[1])

                (tn, tnl1), stats = specimen.collect(tn)
                tnl += tnl1
                res.count(stats[0], stats[1])

            (tp1, fn), stats = self.SSM.test_tb(fn)
            tp += tp1
            res.count(stats[0], stats[1])

            (fp1, tn), stats = self.SSM.test_nontb(tn)
            fp += fp1
            res.count(stats[0], stats[1])
            res.count('N_Det_SSM', tp1 + fp1)

        if self.Xpert is not None:
            if self.Xpert.Specimen is not None:
                if specimen is not self.Xpert.Specimen:
                    specimen = self.Xpert.Specimen
                    (fn, fnl), stats = specimen.collect(fn + fnl)
                    res.count(stats[0], stats[1])

                    (tn, tnl), stats = specimen.collect(tn + tnl)
                    res.count(stats[0], stats[1])

            else:
                fn, fnl = fn + fnl, 0
                tn, tnl = tn + tnl, 0

            (tp1, fn), stats = self.Xpert.test_tb(fn)
            tp += tp1
            res.count(stats[0], stats[1])

            (fp1, tn), stats = self.Xpert.test_nontb(tn)
            fp += fp1
            res.count(stats[0], stats[1])

            res.count('N_Det_Xpert', tp1 + fp1)

        res.count('N_SampleFailed', fnl)
        res.count('N_SampleFailed', tnl)

        if self.CDx is not None:
            (tp1, fn), stats = self.CDx.test_tb(fn)
            tp += tp1
            res.count(stats[0], stats[1])

            (tp1, fnl), stats = self.CDx.test_tb(fnl)
            tp += tp1
            res.count(stats[0], stats[1])
            res.count('N_PreDxLTFU', fnl)

            (fp1, tn), stats = self.CDx.test_nontb(tn)
            fp += fp1
            res.count(stats[0], stats[1])

            (fp1, tnl), stats = self.CDx.test_nontb(tnl)
            fp += fp1
            res.count(stats[0], stats[1])
            res.count('N_PreDxLTFU', tnl)

        res.push(tp=tp, fn=fn + fnl, fn_pd=fnl, fp=fp, tn=tn + tnl, tn_pd=tnl)
        res.count(f'N_Eval_TB_{self.Key}', n_tb)
        res.count(f'N_Eval_NonTB_{self.Key}', n_nontb)
        res.count(f'N_Det_TB_{self.Key}', tp)
        res.count(f'N_Det_NonTB_{self.Key}', fp)
        return res

    def dx_sto(self, n_tb=0, n_nontb=0):
        res = Results()

        tp, fn, fnl = 0, n_tb, 0
        fp, tn, tnl = 0, n_nontb, 0

        specimen = None
        if self.SSM is not None:
            specimen = self.SSM.Specimen
            if specimen is not None:
                (fn, fnl1), stats = specimen.collect_sto(fn)
                fnl += fnl1
                res.count(stats[0], stats[1])

                (tn, tnl1), stats = specimen.collect_sto(tn)
                tnl += tnl1
                res.count(stats[0], stats[1])

            (tp1, fn), stats = self.SSM.test_tb_sto(fn)
            tp += tp1
            res.count(stats[0], stats[1])

            (fp1, tn), stats = self.SSM.test_nontb_sto(tn)
            fp += fp1
            res.count(stats[0], stats[1])
            res.count('N_Det_SSM', tp1 + fp1)

        if self.Xpert is not None:
            if self.Xpert.Specimen is not None:
                if specimen is not self.Xpert.Specimen:
                    specimen = self.Xpert.Specimen
                    (fn, fnl), stats = specimen.collect_sto(fn + fnl)
                    res.count(stats[0], stats[1])

                    (tn, tnl), stats = specimen.collect_sto(tn + tnl)
                    res.count(stats[0], stats[1])

            else:
                fn, fnl = fn + fnl, 0
                tn, tnl = tn + tnl, 0

            (tp1, fn), stats = self.Xpert.test_tb_sto(fn)
            tp += tp1
            res.count(stats[0], stats[1])

            (fp1, tn), stats = self.Xpert.test_nontb_sto(tn)
            fp += fp1
            res.count(stats[0], stats[1])
            res.count('N_Det_Xpert', tp1 + fp1)

        if self.CDx is not None:
            (tp1, fn), stats = self.CDx.test_tb_sto(fn)
            tp += tp1
            res.count(stats[0], stats[1])

            (tp1, fnl), stats = self.CDx.test_tb_sto(fnl)
            tp += tp1
            res.count(stats[0], stats[1])
            res.count('N_PreDxLTFU', fnl)

            (fp1, tn), stats = self.CDx.test_nontb_sto(tn)
            fp += fp1
            res.count(stats[0], stats[1])

            (fp1, tnl), stats = self.CDx.test_nontb_sto(tnl)
            fp += fp1
            res.count(stats[0], stats[1])
            res.count('N_PreDxLTFU', tnl)

        res.push(tp=tp, fn=fn + fnl, fn_pd=fnl, fp=fp, tn=tn + tnl, tn_pd=tnl)
        res.count(f'N_Eval_TB_{self.Key}', n_tb)
        res.count(f'N_Eval_NonTB_{self.Key}', n_nontb)
        res.count(f'N_Det_TB_{self.Key}', tp)
        res.count(f'N_Det_NonTB_{self.Key}', fp)
        return res


if __name__ == '__main__':

    print('Test')
    t0 = Test('SSM', 0.8, 0.98)

    print('-- TB')
    print(t0.test_tb(1000))
    print(t0.test_tb_sto(1000))

    print('-- NonTB')
    print(t0.test_nontb(1000))
    print(t0.test_nontb_sto(1000))

    sputum = Specimen('Sputum', 0.1)
    ssm = Test('SSM', 0.6, 0.98, sputum)
    xpert = Test('Xpert', 0.8, 0.99, sputum)
    cdx = Test('CDx', 0.7, 0.95)

    alg1 = Algorithm('SSM > Xpert > CDx', ssm=ssm, xpert=xpert, cdx=cdx)
    alg2 = Algorithm('SSM > CDx', ssm=ssm, cdx=cdx)
    alg3 = Algorithm('Xpert > CDx', xpert=xpert, cdx=cdx)
    alg4 = Algorithm('CDx', cdx=cdx)

    res = Results()
    for alg in [alg1, alg2, alg3, alg4]:
        print(' ')
        print(alg.Key)
        res1 = alg.dx(1000, 10000)
        res1.print(True)
        res += res1

        print(' ')
        alg.dx_sto(1000, 10000).print(True)

    print('\n\n')
    res.print(detail=True)
