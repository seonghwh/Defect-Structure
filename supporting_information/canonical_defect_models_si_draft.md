# Supporting Information: Canonical Defect Model Used for BaTiO3, Mn-BaTiO3, and Mn+Y-BaTiO3

## S1. Purpose and scope of the model

Canonical defect models were used to interpret the relative changes in electronic carrier concentration, oxygen vacancy concentration, and dopant charge state among the BaTiO3-based ceramics examined in the main text. The model follows the general BaTiO3 defect-chemical framework developed in prior studies of intrinsic, acceptor-doped, Mn-doped, and Mn+Y co-doped BaTiO3 [1-11]. The model was not intended to provide a complete atomistic description of every defect complex in the experimental ceramics. Instead, it provides a constrained defect-chemical framework for comparing the four experimental compositions under the same thermal and atmospheric conditions:

1. nominally undoped BaTiO3 with A/B = 0.999,
2. 0.5 mol% Mn-doped BaTiO3 with A/B = 0.994,
3. 0.5 mol% Mn-doped BaTiO3 with A/B = 0.999, and
4. 0.5 mol% Mn+Y co-doped BaTiO3 with A/B = 0.999.

The calculated defect diagrams in the main text were obtained by first equilibrating the cation and oxygen defect populations at 1150 C and then calculating the low-temperature state at 500 C. The cation vacancy concentrations were frozen from the 1150 C equilibrium state, while electronic carriers, oxygen vacancies, and Mn redox states were allowed to re-equilibrate at 500 C as a function of oxygen partial pressure. This treatment approximates the fact that cation diffusion is slow during cooling, whereas oxygen exchange and electronic redistribution remain more accessible over the temperature range relevant to the impedance measurements.

All concentrations are reported in cm^-3. Oxygen partial pressure is expressed in atm. The Boltzmann constant used in the calculations was k = 8.61733e-5 eV K^-1. The B-site density used to convert A/B nonstoichiometry into a fixed A-site vacancy concentration was N_B = 1.55e22 cm^-3.

## S2. Defect species and notation

The following species were included explicitly:

| Symbol in model | Kroger-Vink interpretation | Role in charge neutrality |
|:--|:--|:--|
| n | electron, e' | negative electronic carrier |
| p | hole, h* | positive electronic carrier |
| V_O** | doubly ionized oxygen vacancy | positive ionic defect |
| V_Ba'' | doubly charged barium vacancy | negative cation vacancy |
| V_Ti'''' | quadruply charged titanium vacancy | negative cation vacancy |
| Mn_Ti^x | tetravalent Mn on Ti site, Mn4+ | neutral relative to Ti4+ |
| Mn_Ti' | trivalent Mn on Ti site, Mn3+ | singly negative acceptor state |
| Mn_Ti'' | divalent Mn on Ti site, Mn2+ | doubly negative acceptor state |
| Y_Ba* | Y donor on Ba site | singly positive donor state |
| A' | background acceptor in nominally undoped BaTiO3 | singly negative acceptor |

For nominally undoped BaTiO3, a small background acceptor concentration of 100 ppm relative to N_B was included, giving [A'] = 1.55e18 cm^-3. This term accounts for residual acceptor impurities in nominally undoped BaTiO3 and avoids treating the commercial starting powder as an ideal intrinsic crystal, consistent with prior discussions that nominally undoped BaTiO3 can be strongly influenced by extrinsic impurities and nonstoichiometry [1,2,5].

## S3. Defect equilibria

### S3.1 Intrinsic electronic disorder

The electron-hole equilibrium was written as

```text
nil <=> e' + h*
```

with

```text
K_i = n p
p = K_i / n
K_i = 8.55e44 exp(-2.91/kT).
```

The same K_i was used for all compositions.

### S3.2 Oxygen reduction equilibrium

The oxygen reduction reaction was written as

```text
O_O^x <=> V_O** + 2e' + 1/2 O2(g)
```

with the concentration-based mass-action expression

```text
K_R = [V_O**] n^2 pO2^(1/2)
[V_O**] = K_R / (n^2 pO2^(1/2)).
```

Two effective values of K_R were used:

| Composition class | K_R expression |
|:--|:--|
| nominally undoped BaTiO3 | K_R = 2.56e71 exp(-6.10/kT) |
| Mn-containing BaTiO3 and Mn+Y-BaTiO3 | K_R = 1.06e71 exp(-5.69/kT) |

The reason for using two K_R values is discussed in Section S5.

### S3.3 Schottky disorder and cation nonstoichiometry

For stoichiometric A/B = 1 compositions, Schottky disorder was represented as

```text
nil <=> V_Ba'' + V_Ti'''' + 3V_O**
```

with

```text
K_S = [V_Ba''][V_Ti''''][V_O**]^3
K_S = 3.4e105 exp(-2.795/kT).
```

For A/B = 1, the model assumed [V_Ba''] = [V_Ti'''']. Combining the Schottky equilibrium with the oxygen reduction expression gives

```text
[V_Ba''] = [V_Ti''''] = sqrt(K_S / K_R^3) n^3 pO2^(3/4).
```

For Ti-rich compositions with A/B < 1, the A-site vacancy concentration was fixed by the nominal cation ratio:

```text
[V_Ba''] = (1 - A/B) N_B
[V_Ti''''] = 0.
```

Thus A/B = 0.999 gives [V_Ba''] = 1.55e19 cm^-3, and A/B = 0.994 gives [V_Ba''] = 9.30e19 cm^-3.

### S3.4 Mn redox equilibria

Mn was assumed to occupy the Ti site and to access Mn4+, Mn3+, and Mn2+ charge states. The two electron exchange reactions were

```text
Mn_Ti' <=> Mn_Ti^x + e'
Mn_Ti'' <=> Mn_Ti' + e'
```

with

```text
K_Mn43 = n[Mn_Ti^x]/[Mn_Ti']
K_Mn43 = 3.2e22 exp(-1.84/kT)

K_Mn32 = n[Mn_Ti']/[Mn_Ti'']
K_Mn32 = 0.8e22 exp(-1.25/kT).
```

The total Mn concentration was fixed at 7.7e19 cm^-3 for the Mn-doped and Mn+Y co-doped compositions. The Mn site balance was

```text
Mn_total = [Mn_Ti^x] + [Mn_Ti'] + [Mn_Ti''].
```

For a given electron concentration n, the Mn partitioning was calculated as

```text
[Mn_Ti']  = Mn_total / (K_Mn43/n + 1 + n/K_Mn32)
[Mn_Ti''] = [Mn_Ti'] n / K_Mn32
[Mn_Ti^x] = Mn_total - [Mn_Ti'] - [Mn_Ti''].
```

### S3.5 Y donor concentration

For the Mn+Y co-doped composition, Y was treated as a fixed donor on the Ba site:

```text
[Y_Ba*] = Y_total = 7.7e19 cm^-3.
```

This treatment follows the donor-acceptor compensation picture for Mn+Y co-doped BaTiO3, in which Y_Ba* provides a donor reservoir while the multivalent Mn acceptor changes charge state with oxygen activity and Fermi level position [10,11].

## S4. Charge neutrality and numerical solution

The charge neutrality condition was

```text
n + 2[V_Ba''] + 4[V_Ti''''] + [Mn_Ti'] + 2[Mn_Ti''] + z_acc[Acc]
    = p + 2[V_O**] + [Y_Ba*].
```

Here z_acc = 1 for the singly charged background acceptor A'. For each value of pO2, the equation was solved numerically for n in log10(n) space. Once n was obtained, p, [V_O**], [V_Ba''], [V_Ti''''], and the Mn charge-state populations were calculated from the mass-action and site-balance expressions above.

The pO2 grid used for the defect diagrams ranged from 1e-20 to 1e5 atm. The reported air values were calculated at pO2 = 0.21 atm, not by interpolation from the grid.

## S5. Rationale for using two effective K_R models

For an ideal thermodynamic reaction written strictly in terms of activities, the standard equilibrium constant for oxygen reduction should be unique at a given temperature. The present model, however, uses a simplified concentration-based mass-action expression,

```text
K_R = [V_O**] n^2 pO2^(1/2),
```

and does not explicitly include activity coefficients, dopant-vacancy association, local charge compensation, space-charge effects, or polaronic carrier trapping. Therefore the fitted K_R in this model is an apparent or effective reduction constant. It contains not only the standard oxygen-reduction thermodynamics but also the nonideal defect interactions and activity corrections that are omitted from the canonical concentration-only expression.

This distinction is important for Mn-containing BaTiO3. Mn is not a fixed-valence acceptor. It can exist as Mn4+, Mn3+, and Mn2+ on the Ti site, and the relative populations depend on pO2 and temperature [5-9]. Hagemann and Hennings showed by reversible weight-change measurements that acceptor-doped BaTiO3 has oxygen vacancy concentrations and acceptor valence states that depend strongly on annealing pO2 and temperature, and that the charge deficiency of acceptors is compensated by doubly ionized oxygen vacancies [5]. In particular, their Mn-doped samples showed changes in Mn valence state accompanying oxygen vacancy formation. Thus the oxygen reduction process in Mn-doped BaTiO3 is coupled to dopant redox chemistry rather than being identical to the nominally undoped case.

Later work reached the same physical conclusion from transport, nonstoichiometry, and spectroscopic perspectives. Song, Yoo, and Kim analyzed Mn-doped BaTiO3 transport in terms of Mn ionization equilibria [6], while Lee, Yoo, and Becker reported that the oxygen nonstoichiometry of Mn-doped BaTiO3 is larger than that of undoped BaTiO3 and that Mn changes valence from Mn4+ to Mn3+ to Mn2+ with decreasing pO2 [7]. Yoon, Randall, and Hur further showed that variable-valence Mn-doped BaTiO3 differs strongly from fixed-valence Mg-doped BaTiO3 during cooling and reoxidation: Mn-doped BaTiO3 has a much lower oxygen vacancy concentration than Mg-doped BaTiO3 under comparable low-temperature oxidizing conditions, together with stronger electron/hole trapping [8]. Chikada et al. directly supported this interpretation by combining first-principles calculations and ESR measurements, showing that oxygen vacancies are stabilized next to Mn and that electrons generated by oxygen vacancy formation are trapped in Mn 3d states [9]. These effects are precisely the types of nonideal interactions that are folded into an effective K_R when they are not modeled explicitly.

The use of two effective K_R values is therefore a controlled approximation. The model does not independently tune K_R for every sample. Instead, it uses one effective reduction constant for the nominally undoped BaTiO3 baseline and one Mn-containing effective reduction constant for both Mn-doped and Mn+Y co-doped BaTiO3. The A/B ratio, background acceptor concentration, Mn redox equilibria, and Y donor concentration are then treated explicitly. This keeps the model constrained while allowing it to capture the experimentally required distinction between the nominally undoped defect chemistry and the Mn-bearing defect chemistry.

The Mn+Y co-doped composition uses the same Mn-containing K_R as the Mn-only compositions because the additional Y effect is included explicitly through [Y_Ba*] in the charge neutrality equation. This choice separates the Mn-related effective oxygen-reduction behavior from the donor-acceptor compensation effect of Y. The resulting calculation reproduces the expected suppression of both hole and oxygen vacancy concentrations in the Mn+Y sample, consistent with the Fermi-level-pinning picture reported for BaTiO3 co-doped with a variable-valence acceptor and a fixed-valence donor [10,11].

## S6. Equilibrium constants used in the final model

| Constant | Expression | Used for | Source / basis |
|:--|:--|:--|:--|
| K_i | 8.55e44 exp(-2.91/kT) | all compositions | canonical BaTiO3 electronic disorder framework [1,2] |
| K_S | 3.4e105 exp(-2.795/kT) | all compositions | canonical BaTiO3 Schottky disorder framework [1,2] |
| K_R, undoped effective | 2.56e71 exp(-6.10/kT) | nominally undoped BaTiO3 | undoped baseline reduction behavior used to match the high-conductivity BaTiO3 baseline; related to BaTiO3 reduction-constant values used in prior defect models [3,4,8] |
| K_R, Mn effective | 1.06e71 exp(-5.69/kT) | Mn-doped and Mn+Y co-doped BaTiO3 | Mn-bearing effective reduction behavior, consistent with Mn redox / vacancy interaction literature [5-9] |
| K_Mn43 | 3.2e22 exp(-1.84/kT) | Mn-containing compositions | Mn4+/Mn3+ redox equilibrium [5-9] |
| K_Mn32 | 0.8e22 exp(-1.25/kT) | Mn-containing compositions | Mn3+/Mn2+ redox equilibrium [5-9] |

## S7. Calculated air values used in the manuscript

The following values were calculated at T = 500 C and pO2 = 0.21 atm after freezing the cation vacancy concentration from the 1150 C equilibrium state.

| Composition | n | p | V_O** | V_Ba'' | V_Ti'''' | Mn_Ti^x | Mn_Ti' | Mn_Ti'' | Y_Ba* | A' |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Undoped BaTiO3, A/B = 0.999 | 4.188e6 | 2.175e19 | 5.401e18 | 1.550e19 | 0 | 0 | 0 | 0 | 0 | 1.550e18 |
| 0.5 mol% Mn-BaTiO3, A/B = 0.994 | 1.434e7 | 6.351e18 | 8.984e19 | 9.300e19 | 0 | 7.697e19 | 3.421e16 | 8.657e9 | 0 | 0 |
| 0.5 mol% Mn-BaTiO3, A/B = 0.999 | 3.598e7 | 2.532e18 | 1.428e19 | 1.550e19 | 0 | 7.691e19 | 8.576e16 | 5.444e10 | 0 | 0 |
| 0.5 mol% Mn+Y-BaTiO3, A/B = 0.999 | 4.775e10 | 1.908e15 | 8.106e12 | 1.550e19 | 0 | 3.104e19 | 4.592e19 | 3.869e16 | 7.700e19 | 0 |

## S8. Code availability and reproducibility

The model equations were implemented in Python. The final model is contained in `src/batio3_defects/solver.py`. The defect diagrams corresponding to the manuscript Figure 6 are generated by

```text
python scripts/plot_figure6_defect_concentrations.py
```

This script writes the plotted data to `results/figure6_defect_concentrations_data.csv` and the figure files to `figures/figure6_defect_concentrations.png`, `.pdf`, and `.svg`.

## References

1. J. Nowotny and M. Rekas, "Defect chemistry of BaTiO3," Solid State Ionics 49, 135-154 (1991). https://doi.org/10.1016/0167-2738(91)90079-Q

2. J. Nowotny and M. Rekas, "Defect structure, electrical properties and transport in barium titanate. VI. General defect model," Ceramics International 20, 257-263 (1994).

3. J. Daniels and K. H. Hardtl, "Electrical conductivity at high temperatures of donor-doped barium titanate ceramics. I," Philips Research Reports 31, 489-504 (1976).

4. J. Daniels, "Electrical conductivity at high temperatures of donor-doped barium titanate ceramics. II," Philips Research Reports 31, 505-515 (1976).

5. H. J. Hagemann and D. Hennings, "Reversible weight change of acceptor-doped BaTiO3," Journal of the American Ceramic Society 64, 590-594 (1981). https://doi.org/10.1111/j.1151-2916.1981.tb10223.x

6. C. R. Song, H. I. Yoo, and J. Y. Kim, "Mn-doped BaTiO3: Electrical transport properties in equilibrium state," Journal of Electroceramics 1, 27-39 (1997). https://doi.org/10.1023/A:1009994230779

7. D. K. Lee, H. I. Yoo, and K. D. Becker, "Nonstoichiometry and defect structure of Mn-doped BaTiO3," Solid State Ionics 154-155, 189-193 (2002).

8. S. H. Yoon, C. A. Randall, and K. H. Hur, "Difference between resistance degradation of fixed valence acceptor (Mg) and variable valence acceptor (Mn)-doped BaTiO3 ceramics," Journal of Applied Physics 108, 064101 (2010). https://doi.org/10.1063/1.3480992

9. S. Chikada, T. Kubota, A. Honda, S. Higai, Y. Motoyoshi, N. Wada, and K. Shiratsuyu, "Interactions between Mn dopant and oxygen vacancy for insulation performance of BaTiO3," Journal of Applied Physics 120, 142122 (2016). https://doi.org/10.1063/1.4963381

10. Y. Y. Yeoh, H. Jang, and H. I. Yoo, "Defect structure and Fermi-level pinning of BaTiO3 co-doped with a variable-valence acceptor (Mn) and a fixed-valence donor (Y)," Physical Chemistry Chemical Physics 14, 1642-1648 (2012). https://doi.org/10.1039/C2CP22711H

11. C. E. Lee, S. H. Kang, D. S. Sinn, and H. I. Yoo, "Co-doping effect of Mn and Y on charge and mass transport properties of BaTiO3," Journal of Electroceramics 13, 785-791 (2004). https://doi.org/10.1007/s10832-004-5193-9
