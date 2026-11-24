# PV-LAND offline

Modèle de surface PV-LAND (parc photovoltaïque au sol) forcé par des mesures
météorologiques. La physique est **celle d'ARPS** (`arps5.3.4_PVLAND`, option
`sfcphy = 5`) : ARPS est la référence, l'offline n'est qu'un pilote qui appelle
exactement les mêmes sources.

## Organisation

| Répertoire | Contenu | Modifiable ici ? |
|---|---|---|
| `core/` | physique PV-LAND : copies de `src/arps/pvland*.f90` et `module_pvland.f90` d'ARPS | **non** |
| `arps_lib/` | routines ARPS utilisées par PV-LAND (`soilebm3d`, `thermolib3d`, `sfcflx`/`dragcoef`… extraits de `sfcphy3d.f90`, `tridiag2` de `tmix3d.f90`) et fichiers `.inc` | **non** |
| `offline/pvland_sfcphy_offline.f90` | un pas de temps PV-LAND : **généré** à partir du bloc `sfcphy == 5` de `sfcphy3d.f90` d'ARPS | **non** |
| `offline/pvland_offline_main.f90` | lecture du CSV, état atmosphérique en colonne, boucle temporelle, sorties | oui |
| `offline/pvland_park_map.f90` | emplacement du parc (tout le domaine en offline) | oui |
| `scripts/` | `sync_from_arps.py` (copie/génération depuis ARPS ou vérification) et le gabarit `pvland_sfcphy_offline.tmpl.f90` | oui |
| `legacy/` | anciennes versions non compilées | — |

Le pilote généré ne diffère du bloc ARPS que par les lignes marquées `[OFFLINE]` :
l'offline utilise le **rayonnement diffus mesuré** (au lieu de `split_sw`), les
**angles λ, φ mesurés** (au lieu de `compute_lambda_phi`) et la **longueur
d'Obukhov mesurée** (au lieu de celle calculée à partir des flux). Les lignes
ARPS d'origine restent visibles en commentaire juste au-dessus.

## Compilation et exécution

```bash
make              # compile pv_land.exe (ifx -O2)
make run          # forçage complet (été 2025, 2303 pas de 30 min, ~2 min)
make run-short    # 5 premiers jours (recorded_data_short_5days.csv, ~15 s)
./pv_land.exe mon_forcage.csv
```

Sorties : `pv_land_output.csv` (flux modélisés et mesurés, températures) et
`pv_diag_energy_balance.csv` (bilan d'énergie des panneaux).

## Modifier PV-LAND

1. Faire la modification dans ARPS (`arps5.3.4_PVLAND/src/arps/`).
2. `make sync-arps` : recopie `core/`, `arps_lib/` et régénère le pilote.
3. `make && make run-short` pour tester.

On peut aussi modifier un fichier de `core/` ici pour tester, puis le recopier
dans ARPS ; `make check-arps` signale tout fichier qui diffère d'ARPS.

Les paramètres du parc (géométrie, propriétés optiques) sont dans
`INIT_PVLAND` (`core/pvland.f90`), comme dans ARPS. Les paramètres du site
(sol, végétation, état initial) sont en tête de `offline/pvland_offline_main.f90`.
