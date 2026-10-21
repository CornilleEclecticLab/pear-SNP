#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s12.x1.downsample_AD_count.py
# @Description:
#    In silico downsampling of allelic read counts.
#    It would be launched by s12.downsample_AD_count_01_count.sh

import argparse
import sys
import zlib
import numpy as np

ERR = 0.01
LOG_OK, LOG_ERR, LOG_HALF = np.log10(1 - ERR), np.log10(ERR), np.log10(0.5)
CHUNK = 20000


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--chrom', required=True)
    p.add_argument('--samples', required=True, help='bcftools query -l output (column order)')
    p.add_argument('--depth', required=True, help='Sample_ID<TAB>mapped depth')
    p.add_argument('--sift', required=True, help='merged.SIFTannotations.xls.cleaned')
    p.add_argument('--target', type=float, default=10.0)
    p.add_argument('--reps', type=int, default=3)
    p.add_argument('--seed', type=int, default=20260923)
    p.add_argument('--out', required=True)
    return p.parse_args()


def load_sift(file_path, chrom):
    ann = {}
    with open(file_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if parts[0] != chrom:
                continue
            variant_type, sift_score, prediction = parts[8], parts[12], parts[16]
            if sift_score == 'NA':
                continue
            if 'Low confidence' in prediction or 'WARNING' in prediction:
                raise ValueError(f'Low confidence SIFT prediction at {chrom}:{parts[1]}')
            if variant_type == 'SYNONYMOUS':
                cls = 'SYN'
            elif variant_type == 'NONSYNONYMOUS':
                cls = prediction
            else:
                continue
            ann[int(parts[1])] = (parts[2], parts[3], cls)
    return ann


def gt_code(gt):
    # -1 missing, 0 hom-ref, 1 het, 2 hom-alt
    alleles = gt.replace('|', '/').split('/')
    n1 = alleles.count('1')
    n0 = alleles.count('0')
    if n0 + n1 == 0:
        return -1
    if n0 + n1 == 1:
        raise ValueError(f'Invalid genotype with missing data: {gt}')
    return n1


def call(r, a, d):
    # r, a, d: int arrays (sites x samples). Return genotype codes with filters applied.
    l0 = r * LOG_OK + a * LOG_ERR
    l1 = (r + a) * LOG_HALF
    l2 = r * LOG_ERR + a * LOG_OK
    pl = -10 * np.stack([l0, l1, l2], axis=-1)
    pl -= pl.min(axis=-1, keepdims=True)
    gt = pl.argmin(axis=-1)
    gq = np.minimum(np.sort(pl, axis=-1)[..., 1], 99)
    return np.where((d < 3) | (d > 100) | (gq < 20), -1, gt)


def main():
    args = parse_args()

    with open(args.samples) as f:
        samples = [s.strip() for s in f if s.strip()]
    depth = {}
    with open(args.depth) as f:
        for line in f:
            sid, dp = line.split()[0:2]
            depth[sid] = float(dp)
    prob = np.array([min(1.0, args.target / depth[s]) for s in samples])

    ann = load_sift(args.sift, args.chrom)
    sys.stderr.write(f'{args.chrom}: {len(ann)} SIFT4G sites\n')

    tag = f'ds{args.target:g}x'
    arms = ['orig', 'recall_full'] + [f'{tag}_rep{k}' for k in range(1, args.reps + 1)]
    rngs = [np.random.default_rng([args.seed, k, zlib.crc32(args.chrom.encode())])
            for k in range(1, args.reps + 1)]
    ns = len(samples)
    # counts[arm] = [del_all, del_ho, del_he, syn_site] per sample
    counts = {arm: np.zeros((4, ns), dtype=np.int64) for arm in arms}
    n_sites = {'SYN': 0, 'DELETERIOUS': 0, 'TOLERATED': 0, 'mismatch': 0}

    def flush(cls_list, G, R, A, D):
        if not cls_list:
            return
        G, R, A, D = (np.array(x, dtype=np.int64) for x in (G, R, A, D))
        is_del = np.array([c == 'DELETERIOUS' for c in cls_list])[:, None]
        is_syn = np.array([c == 'SYN' for c in cls_list])[:, None]
        valid = G >= 0   # genotypes missing in the original VCF stay missing

        genos = {'orig': G, 'recall_full': np.where(valid, call(R, A, D), -1)}
        for k, rng in enumerate(rngs, start=1):
            r2 = rng.binomial(R, prob)
            a2 = rng.binomial(A, prob)
            d2 = r2 + a2 + rng.binomial(np.maximum(D - R - A, 0), prob)
            genos[f'{tag}_rep{k}'] = np.where(valid, call(r2, a2, d2), -1)

        for arm, g in genos.items():
            c = counts[arm]
            c[0] += (is_del * ((g == 2) * 2 + (g == 1))).sum(axis=0)
            c[1] += (is_del * (g == 2) * 2).sum(axis=0)
            c[2] += (is_del * (g == 1)).sum(axis=0)
            c[3] += (is_syn * (g >= 1)).sum(axis=0)

    cls_list, G, R, A, D = [], [], [], [], []
    seen = set()
    for line in sys.stdin:
        parts = line.rstrip('\n').split('\t')
        pos = int(parts[1])
        if pos in seen:       # -R with regions may repeat a record
            continue
        seen.add(pos)
        if pos not in ann:
            continue
        ref, alt, cls = ann[pos]
        if parts[2] != ref or parts[3] != alt:
            n_sites['mismatch'] += 1
            continue
        n_sites[cls] += 1
        g_row, r_row, a_row, d_row = [], [], [], []
        for field in parts[4:]:
            gt, ad, dp = field.split(':')
            code = gt_code(gt)
            ads = ad.split(',')
            if code >= 0 and ads[0] != '.' and len(ads) == 2:
                r_row.append(int(ads[0]))
                a_row.append(int(ads[1]))
                d_row.append(int(dp) if dp != '.' else int(ads[0]) + int(ads[1]))
            else:
                r_row.append(0); a_row.append(0); d_row.append(0)
            g_row.append(code)
        cls_list.append(cls); G.append(g_row); R.append(r_row); A.append(a_row); D.append(d_row)
        if len(cls_list) >= CHUNK:
            flush(cls_list, G, R, A, D)
            cls_list, G, R, A, D = [], [], [], [], []
    flush(cls_list, G, R, A, D)

    sys.stderr.write(f'{args.chrom}: sites used {n_sites}\n')
    with open(args.out, 'w') as fo:
        fo.write('Arm\tChrom\tSample_ID\tDel_all\tDel_ho\tDel_he\tSyn_site\n')
        for arm in arms:
            c = counts[arm]
            for i, s in enumerate(samples):
                fo.write(f'{arm}\t{args.chrom}\t{s}\t{c[0, i]}\t{c[1, i]}\t{c[2, i]}\t{c[3, i]}\n')


if __name__ == '__main__':
    main()
