# PV-LAND offline — version 2

Modèle de surface PV-LAND (parc photovoltaïque au sol) forcé par des mesures
météorologiques. La physique est **celle d'ARPS** (`arps5.3.4_PVLAND_v2`,
branche `v2` de arps_pvland, option `sfcphy = 5`) : ARPS est la référence,
l'offline n'est qu'un pilote qui appelle exactement les mêmes sources.

La version 2 reprend la physique de `PV_LAND_final` (2025) dans la structure
de la version 1 (branche `main`, `PVLAND_24july2026`) :

- longueur d'Obukhov calculée par le modèle (`pvland_obukhov.f90` : flux de
  flottabilité H + LE, relaxation sur `pv_tau_L`, |L| ≥ `pv_L_floor`, L = 0
  neutre) ; le forçage ne contient plus de L mesurée ;
- humidité de l'air sous les panneaux `q_can` (nœud sans stockage,
  conductance `gcan`), F3 d'Alfieri et al. (2008) sur le déficit de saturation
  de cet air ;
- sol propre à PV-LAND (`pvland_soil.f90`) : drainage, ruissellement de la
  canopée, seuil `pv_wref`, paramètres de sol et `cgv` dans la namelist ;
- coefficients d'échange `dragcoef_pv` (`pvland_dragcoef.f90`) avec zt/z0
  explicite (`pv_ztz0_nat`, `pv_ztz0_air`) ;
- convection des panneaux h = C·U^p (`pv_C_for_*`, `pv_p_for_*`), convection
  libre ×`pv_fac_h_free`, plus de facteurs multiplicatifs ;
- profil de vent dépendant de `z0_nat` (= `roufns`) et `pv_zd_nat` ;
- facteurs de vue, ombre au sol exacte, fractions éclairées des panneaux en
  tan λ·sin φ (λ = angle zénithal, φ = azimut dans le repère du parc) ;
- paramètres du module PV et émissivités de PV_LAND_final.

Choix gardés d'ARPS : simple précision, H_nat = flux de `soilebm` (`shflx`),
z_ref tiré de `zp`, vent minimal `vsfcmin`. Les routines générales d'ARPS
(`soilebm`, `dragcoef`, `soilcst.inc`) ne sont pas modifiées : `sfcphy = 3/4`
restent inchangés.

Validation (forçage 2025, 2350 pas) : écarts moyens à PV_LAND_final
< 0,1 W/m² sur H, LE, Rn, G et < 0,02 K sur les températures ; RMSE face aux
mesures identiques (H 47,8, LE 55,7, Rn 13,5, G 10,1 W/m²).

## Organisation

| Répertoire | Contenu | Modifiable ici ? |
|---|---|---|
| `core/` | physique PV-LAND : copies de `src/arps/pvland*.f90` et `module_pvland.f90` d'ARPS, dont `pvland_driver.f90` (un pas de temps PV-LAND, `pvland_sfcphy`, appelé aussi par ARPS), `pvland_obukhov.f90`, `pvland_soil.f90`, `pvland_dragcoef.f90` | **non** |
| `arps_lib/` | routines ARPS utilisées par PV-LAND (`soilebm3d`, `thermolib3d`, `sfcflx`/`dragcoef`… extraits de `sfcphy3d.f90`, `tridiag2` de `tmix3d.f90`) et fichiers `.inc` | **non** |
| `offline/pvland_offline_main.f90` | lecture du CSV, état atmosphérique en colonne, boucle temporelle, sorties | oui |
| `offline/pvland_park_map.f90` | emplacement du parc (tout le domaine en offline) | oui |
| `scripts/` | `sync_from_arps.py` (copie depuis ARPS ou vérification), `compare_outputs.py` (comparaison de deux sorties et RMSE face aux mesures) | oui |
| `analysis/` | scripts Python de PV_LAND_final (statistiques, figures ; mesures dans `measured_data_2025.csv`) | oui |
| `legacy/` | anciennes versions non compilées | — |

ARPS (`sfcphy3d.f90`, `sfcphy = 5`) et l'offline appellent la **même routine**
`pvland_sfcphy`. Seul l'argument `meas_opt` diffère : ARPS passe `0` (part
diffuse par `split_sw`, angles λ, φ par `compute_lambda_phi`) ; l'offline passe
`1` et fournit le **rayonnement diffus** et les **angles λ, φ mesurés**. La
longueur d'Obukhov est calculée par PV-LAND dans les deux cas.

## Compilation et exécution

```bash
make              # compile pv_land.exe (ifx -O2)
make run          # forçage complet (recorded_data_2025.csv, 2351 pas de 30 min, ~10 s)
make run-short    # 5 premiers jours (recorded_data_2025_short_5days.csv, ~1 s)
./pv_land.exe mon_forcage.csv                 # paramètres lus dans pvland.nml
./pv_land.exe mon_forcage.csv mes_params.nml  # autre fichier de paramètres
```

Sorties (mêmes noms et colonnes que PV_LAND_final, lues par `analysis/`) :
`pv_land_output_2025.csv` (flux modélisés et mesurés, températures, L, q_can,
z0/d du parc, u*…) et `pv_diag_energy_balance_2025.csv` (bilan d'énergie des
panneaux). Pour comparer à une autre sortie :
`python3 scripts/compare_outputs.py pv_land_output_2025.csv autre.csv --meas`.

## Modifier PV-LAND

1. Faire la modification dans ARPS (`arps5.3.4_PVLAND_v2/src/arps/`).
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
- `pv_park_opt`, `pv_park_i1..j2` : emplacement du parc sur la grille ARPS.
  En offline, `pv_park_opt = 2` fait tourner la colonne **sans panneaux**
  (sol naturel seul, comme `sfcphy = 4` dans ARPS) : le rayonnement net du sol
  est alors calculé avec `pv_albedo_nat` et `pv_emiss_nat`, et les colonnes
  « panneaux » des sorties n'ont plus de sens ;
- géométrie (`pv_L`, `pv_a`, `pv_beta_deg`, `pv_h_PV`, `pv_theta_0_deg`,
  `pv_az_park`), rugosités, albédos et émissivités, longueurs d'interface,
  paramètres thermiques et électriques du module PV.

Les paramètres du site et du run offline sont dans une seconde namelist,
`&offline_site`, à la fin de `pvland.nml` (non lue par ARPS) : fichier de
forçage, pas de temps (`dt_forcing`, `dtsfc`), sol et végétation (`styp`,
`vtyp`, `lai0`, `roufns0`, `veg0`, mêmes noms que dans `&soil_ebm` d'ARPS) et
état initial du sol. Seule la grille (une colonne, `nx = ny = 2`, `nz = 3`)
reste fixée dans `offline/pvland_offline_main.f90`.
