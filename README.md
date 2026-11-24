# PV-LAND offline

Modèle de surface PV-LAND (parc photovoltaïque au sol) forcé par des mesures
météorologiques. La physique est **celle d'ARPS** (`arps5.3.4_PVLAND`, option
`sfcphy = 5`) : ARPS est la référence, l'offline n'est qu'un pilote qui appelle
exactement les mêmes sources.

## Organisation

| Répertoire | Contenu | Modifiable ici ? |
|---|---|---|
| `core/` | physique PV-LAND : copies de `src/arps/pvland*.f90` et `module_pvland.f90` d'ARPS, dont `pvland_driver.f90` (un pas de temps PV-LAND, `pvland_sfcphy`, appelé aussi par ARPS) | **non** |
| `arps_lib/` | routines ARPS utilisées par PV-LAND (`soilebm3d`, `thermolib3d`, `sfcflx`/`dragcoef`… extraits de `sfcphy3d.f90`, `tridiag2` de `tmix3d.f90`) et fichiers `.inc` | **non** |
| `offline/pvland_offline_main.f90` | lecture du CSV, état atmosphérique en colonne, boucle temporelle, sorties | oui |
| `offline/pvland_park_map.f90` | emplacement du parc (tout le domaine en offline) | oui |
| `scripts/` | `sync_from_arps.py` (copie depuis ARPS ou vérification) | oui |
| `legacy/` | anciennes versions non compilées | — |

ARPS (`sfcphy3d.f90`, `sfcphy = 5`) et l'offline appellent la **même routine**
`pvland_sfcphy`. Seul l'argument `meas_opt` diffère : ARPS passe `0` (part
diffuse par `split_sw`, angles λ, φ par `compute_lambda_phi`, longueur
d'Obukhov calculée à partir des flux du pas précédent) ; l'offline passe `1` et
fournit le **rayonnement diffus**, les **angles λ, φ** et la **longueur
d'Obukhov mesurés**.

## Compilation et exécution

```bash
make              # compile pv_land.exe (ifx -O2)
make run          # forçage complet (été 2025, 2303 pas de 30 min, ~2 min)
make run-short    # 5 premiers jours (recorded_data_short_5days.csv, ~15 s)
./pv_land.exe mon_forcage.csv                 # paramètres lus dans pvland.nml
./pv_land.exe mon_forcage.csv mes_params.nml  # autre fichier de paramètres
```

Sorties : `pv_land_output.csv` (flux modélisés et mesurés, températures) et
`pv_diag_energy_balance.csv` (bilan d'énergie des panneaux).

## Modifier PV-LAND

1. Faire la modification dans ARPS (`arps5.3.4_PVLAND/src/arps/`).
2. `make sync-arps` : recopie `core/` et `arps_lib/`.
3. `make && make run-short` pour tester.

On peut aussi modifier un fichier de `core/` ici pour tester, puis le recopier
dans ARPS ; `make check-arps` signale tout fichier qui diffère d'ARPS.

## Paramètres : namelist `&pvland`

Les paramètres du parc sont lus dans la namelist `&pvland` : fichier
`pvland.nml` pour l'offline, bloc `&pvland` de `arps.input` (après `&soil_ebm`)
pour ARPS. C'est la même namelist dans les deux modèles (module `pvland_nml`
de `core/module_pvland.f90`, qui donne aussi les valeurs par défaut).
Chaque paramètre y est décrit sur sa ligne :

- `vf_update_opt` : quand recalculer tout ce qui dépend de l'inclinaison β
  (facteurs de vue, F_sky, factorisations des matrices de radiosité, h_clear,
  z_can, z_PVa…) : `0` = une seule fois à l'initialisation (panneaux fixes),
  `1` = à chaque pas de temps (panneaux mobiles ; la loi de suivi β(t) reste
  à écrire) ;
- `pv_park_opt`, `pv_park_i1..j2` : emplacement du parc sur la grille ARPS
  (non utilisés en offline) ;
- géométrie (`pv_L`, `pv_a`, `pv_beta_deg`, `pv_h_PV`, `pv_theta_0_deg`,
  `pv_az_park`), rugosités, albédos et émissivités, longueurs d'interface,
  paramètres thermiques et électriques du module PV.

Les paramètres du site (sol, végétation, état initial, `dtsfc`) restent en tête
de `offline/pvland_offline_main.f90`.
