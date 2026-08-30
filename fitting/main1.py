# main.py is run by parallel.py

from pathlib import Path
import glob
import numpy as np
import blast

conf = blast.settings


# ----- default worker_name and map function -----

def get_worker_name():
    from datetime import datetime
    date = datetime.now().strftime("%Y%b%d")
    return 'serial_%s' % date


parallel_map = map
blast.runtime['get_worker_name'] = get_worker_name
model_dict = blast.json.read(conf["model"])
blast.runtime["FF"] = blast.model.load_from_dict(model_dict)

# ----- create ram based tmp -----
if 'ramtmp' not in blast.runtime:
    import tempfile
    tmp = tempfile.mkdtemp(prefix="blastff_", dir='/dev/shm')
    if not Path(tmp).exists():
        blast.mpy.ERROR(f"cannot create tmp directory,\n{tmp}")
    else:
        blast.logger.info(f"ram based tmp directory at {tmp}")
    tmp_dir = blast.conf.get('lmp', {}).get('tmp_dir', 'tmp')
    try:
        blast.cmd.ln(tmp, tmp_dir)
    except:
        blast.mpy.ERROR(f"cannot symlink {tmp} to '{tmp_dir}'")
    blast.runtime['ramtmp'] = tmp


# ----- define objective -----
def objective(params):

    runtime = blast.runtime
    wkr = runtime['get_worker_name']
    worker = runtime['worker'][wkr()]

    FF = runtime['FF']

    if 'HO' not in worker:
        HO = worker['HO'] = blast.data.load('ho').setup()
        report_dir = conf.get('data', {}).get('report_dir', 'reports')
        HO.reportfile = f"{report_dir}/ho.report"
        HO.write(str(FF), prefix='# ')
    else:
        HO = worker['HO']

    if 'GEO' not in worker:
        GEO = worker['GEO'] = blast.data.load('structure').catalog()
        GEO.path = blast.conf['data']['path']
    else:
        GEO = worker['GEO']

    tmp_dir = blast.conf.get('lmp', {}).get('tmp_dir', 'tmp')
    wd = worker['wd'] = blast.cwd / f"{tmp_dir}/{wkr()}"
    blast.cmd.mkdir(wd, reset=True)

    HO.reset()
    GEO.reset()

    if FF.get(params) is None:
        return HO.finalobj()
    FF.write(wd / 'trial.ff', verbose=False)
    param_str = ' '.join([str(s) for s in params]) + ' DUMP'

    earlyexit = bool(blast.conf.get('data', {}).get('earlyexit', True))
    if not earlyexit:
        failed = []

    polymorphs_energy_order = sorted([
        Path(fname).name for fname in
        glob.glob(GEO.path + '/*.data')
    ], key=lambda s: int(Path(s).name.split('_')[0]))

    # ----------------------- workflow --------------------

    # PROPERTY CALCULATION: Lattice Parameters
    prop = 'lattice'
    w_stage = np.ones(len(polymorphs_energy_order))
    kwargs = {
        "in|onoff_log": '#',
        "in|replicate": '6 6 1',
        "in|reset": True,
        "in|onoff_boxrelax": '',
        "in|P_couple": 'none',
    }

    for structure in polymorphs_energy_order:
        geo = GEO[structure]

        _t = geo.compute(prop, suffix='_t', verbose=False)
        if _t is SystemExit:
            return HO.finalobj('runtime error')

        _p = geo.compute(prop, suffix='_p', **kwargs, verbose=False)
        if _p is SystemExit:
            return HO.finalobj('runtime error')

        _s = geo.score(prop, tally=True, w_stage=w_stage)
        if _s is SystemExit:
            return HO.finalobj('runtime error')

        _failed = geo.checkpoint(prop, **{"conditions": [
            "values.maxAE% <= 3",
        ]})

        if earlyexit and _failed:
            return HO.finalobj(
                f"{param_str} {structure.split('.')[0]}/"
                f"{polymorphs_energy_order[-1].split('.')[0]} {prop} maxAE > 3%"
            )
        elif not earlyexit:
            failed.append(_failed)

    # PROPERTY CALCULATION: Cohesive Energy and Ordering
    prop = 'ce'
    kwargs = {
        "in|onoff_log": '#',
        "in|reset": True,
        "in|replicate": '6 6 1',
        "in|onoff_boxrelax": '',
        "in|P_couple": 'none',
        "in|onoff_min": '',
    }

    w_stage = np.ones(len(polymorphs_energy_order))
    group = blast.data.load('geo').group()
    for structure in polymorphs_energy_order:
        geo = GEO[structure]

        _t = geo.compute(prop, suffix='_t', verbose=False)
        if _t is SystemExit:
            return HO.finalobj('runtime error')

        for s in _t:
            group._t[s].append(_t[s])

        _p = geo.compute(prop, suffix='_p', **kwargs, verbose=False)
        if _p is SystemExit:
            return HO.finalobj('runtime error')

        for s in _p:
            group._p[s].append(_p[s])

    _s = group.score(
        prop, tally=True, w_stage=w_stage,
        w_typ={'ordering': 0.5, 'values': 1}, typ=['ordering', 'values']
    )
    _failed = group.checkpoint(prop, **{"conditions": [
        "ordering.Norderoff==0",
        "values.maxAE%<=3",
    ]})

    if _s is SystemExit:
        return HO.finalobj('runtime error')

    if earlyexit and _failed:
        return HO.finalobj(
            f"{param_str} {structure.split('.')[0]}/"
            f"{polymorphs_energy_order[-1].split('.')[0]} {prop} maxAE% > 3%"
        )
    elif not earlyexit:
        failed.append(_failed)

    # PROPERTY CALCULATION: Equation of State
    prop = 'eos'

    for structure in polymorphs_energy_order:
        geo = GEO[structure]

        _t = geo.compute(prop, suffix='_t', verbose=False)
        if _t is SystemExit:
            return HO.finalobj('runtime error')

        _p = geo.compute(prop, suffix='_p', strain_basis="a b", verbose=False)
        if _p is SystemExit:
            return HO.finalobj('runtime error')

        _s = geo.score(prop, tally=True)
        if _s is SystemExit:
            return HO.finalobj('runtime error')

        _failed = geo.checkpoint(prop, **{"conditions": [
            "shape.MAE% <= 30",
        ]})

        if earlyexit and _failed:
            return HO.finalobj(
                f"{param_str} {structure.split('.')[0]}/"
                f"{polymorphs_energy_order[-1].split('.')[0]} {prop} failed"
            )
        elif not earlyexit:
            failed.append(_failed)

    # PROPERTY CALCULATION: Elastic Constants
    prop = 'elastic'

    for structure in polymorphs_energy_order:
        geo = GEO[structure]
        geo_prop = geo.meta['prop']
        if 'elastic' not in geo_prop:
            continue

        _t = geo.compute(prop, suffix='_t', verbose=False)
        if _t is SystemExit:
            return HO.finalobj('runtime error')

        _p = geo.compute(prop, suffix='_p', up=1e-6, verbose=False)
        if _p is SystemExit:
            return HO.finalobj('runtime error')

        _s = geo.score(
            prop, tally=True,
            w=np.array([1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        )
        if _s is SystemExit:
            return HO.finalobj('runtime error')

        _failed = geo.checkpoint(prop, **{"conditions": [
            "values.MAE% <= 30"
        ]})

        if earlyexit and _failed:
            return HO.finalobj(
                f"{param_str} {structure.split('.')[0]}/"
                f"{polymorphs_energy_order[-1].split('.')[0]} {prop} MAE% > 30%"
            )
        elif not earlyexit:
            failed.append(_failed)

    # PROPERTY CALCULATION: Phonon Dispersion
    prop = 'phonon'

    structures_with_phonon = [
        s for s in polymorphs_energy_order
        if 'phonon' in GEO[s].meta['prop']
    ]

    w_stage = np.ones(len(structures_with_phonon))
    for structure in structures_with_phonon:
        geo = GEO[structure]

        _t = geo.compute(
            'phonon', suffix='_t',
            png='phonon_t.png',
            txt='phonon_t.txt',
        )
        if _t is SystemExit:
            return HO.finalobj('runtime error')

        _p = geo.compute(
            'phonon', suffix='_p',
            verbose=False,
            png='phonon_p.png',
            txt='phonon_p.txt',
            savedat='phonon_p',
            savedos=True,
            saveCv='Cv_p.txt',
        )
        if _p is SystemExit:
            return HO.finalobj('runtime error')

        _s = geo.score(prop, {**geo._t, **geo._p}, w_stage=w_stage, w_aco=1, w_opt=0)
        if _s is SystemExit:
            return HO.finalobj('runtime error')

        phonon = blast.data.load('phonon')
        phonon.plot({**_t, **_p}, f'plot_{structure}.png')

        _failed = geo.checkpoint('phonon', conditions=[
            "ceil.maxAE% <= 30",
        ])

        if earlyexit and _failed:
            return HO.finalobj(
                f"{param_str} {structure.split('.')[0]}/"
                f"{polymorphs_energy_order[-1].split('.')[0]} phonon failed"
            )
        elif not earlyexit:
            failed.append(_failed)

    return HO.lastfinalobj(f"{param_str} completed all stages")


blast.runtime['objective'] = objective

if __name__ == '__main__':
    import sys
    print(sys.argv[1:])
    objective(sys.argv[1:])
