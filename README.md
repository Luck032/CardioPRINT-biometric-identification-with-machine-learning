# CardioPRINT-biometric-identification-with-machine-learning
This repository contains Python and R programming codes, as well as extracted timestamps for segments that describe emotional states and feature sets for both ECG and ICG recordings that reproduce results for the paper titled "[CardioPRINT: Biometric identification based on the individual characteristics derived from cardiogram](link)" authored by Ilija Tanasković (ORCiD: [0000-0002-6488-4074](https://orcid.org/0000-0002-6488-4074)), Ljiljana B. Lazarević (ORCiD: [0000-0003-1629-3699](https://orcid.org/0000-0003-1629-3699)), Goran Knežević (ORCiD: [0000-0001-8951-3774](https://orcid.org/0000-0001-8951-3774)), Nikola Milosavljević (ORCiD: [0000-0001-5061-149X](https://orcid.org/0000-0001-5061-149X)), Olga Dubljević (ORCiD: [0000-0003-1560-1661](https://orcid.org/0000-0003-1560-1661)), Bojana Bjegojević (ORCiD: [0000-0002-8421-5572](https://orcid.org/0000-0002-8421-5572)) and Nadica Miljković (ORCiD: [0000-0002-3933-6076](https://orcid.org/0000-0002-3933-6076)). Features are calculated  from The dataset was recorded for another study and we share it openly on the Zenodo repository with a Creative Commons Attribution 4.0 International license.

## GitHub Repository Contents
This repository contains Python and R programming codes, as well as extracted timestamps for segments that describe emotional states and feature set for both ECG and ICG recordings that reproduce results for the paper titled "[CardioPRINT: Biometric identification based on the individual characteristics derived from cardiogram](link)".
Also, this repository contains [README.md](link) file with relevant information essential for code reproducibility and [LICENSE](link) file that contains license information that covers shared software codes.

### Code
Shared programs are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. These programs are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with these programs. If not, see [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).

Please, report any bugs to the Authors listed in the Contacts.

The repository contains the following code:

1)	[CardioPrint_feature_extraction_two_segments.R] – R code that implements feature selection on signals containing two segments (emotional phase)
2)	[CardioPrint_feature_extraction_three_segments.R] – R code that implements feature selection on signals containing all three segments
3)	[CardioPrint_all_three_segments_validation.ipynb] – Python code that corresponds to the first step in our Method – The best-performing models and hyperparameter determination (BPM) that implements hyperparameter tuning and the selection of models with high accuracies (>90%) for further analysis
4)	[CardioPrint_all_three_segments_learning.ipynb] – Python code that corresponds to the second step in our Method – Feature Set Selection (FSS) that investigates the best-performing feature set with appropriate statistical tests
5)	[CardioPrint_baseline_learning.ipynb] – Python code that corresponds to the third step in our Method – Testing model robustness to altered emotional state (TMR) that investigates the effect of different emotional states on biometric identification during training model on baseline segment and evaluation on all three segments
6)	[CardioPrint_anger_learning.ipynb] – Python code that corresponds to the third step in our Method – Testing model robustness to altered emotional state (TMR) that investigates the effect of different emotional states on biometric identification during training model on anger segment and evaluation on all three segments

### Data
Data provided in this repository are shared under Attribution 4.0 International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). 

1)	[timestamps_without_neutral.csv] – table with timestamps marking the beginning and the end of each segment in signals containing two segments (emotional phase)
2)	[timestamps_with_neutral.csv] – table with timestamps marking the beginning and the end of each segment in signals containing all three segments (emotional phase)
3)	[Features_without_neutral.csv ] - table containing both ECG and ICG features extracted from signals containing two segments (emotional phase)
4)	[Features.csv] – table containing both ECG and ICG features extracted from signals containing all three segments (emotional phase)
The columns in tables timestamps_without_neutral.csv and  timestamps_with_neutral.csv indicate:
1)	ID – identifier of the individual
2)	T1 – the beginning of the baseline segment
3)	T2 – the end of the baseline segment
4)	T3 – the beginning of the anger segment
5)	T4 – the end of the anger segment
6)	T5 – the beginning of the neutral segment
7)	T6 – the end of the neutral segment
It should be noted that table timestamps_without_neutral.csv does not contain T5 and T6 stamps since the recording protocol did not included neutral segment

## Contacts
Ilija Tanasković ([ilijatanaskovic97@hotmail.com](mailto:ilijatanaskovic97@hotmail.com)) or Nadica Miljković (e-mail: [nadica.miljkovic@etf.bg.ac.rs](mailto:nadica.miljkovic@etf.bg.ac.rs)).

## Funding
NM acknowledges the support from Grant No. 451–03–47/ 2023–01/200103 funded by the Ministry of Science, Technological Development and Innovation of the Republic of Serbia. LBL and GK acknowledge the support from Grant No. 451-03-47/2023-01/200163 funded by the Ministry of Science, Technological Development and Innovation of the Republic of Serbia.

## How to cite this repository?
If you find provided code and signals useful for your own research and teaching class, please cite the following references:

![image](https://github.com/Luck032/CardioPRINT-biometric-identification-with-machine-learning/assets/104569197/f2b01402-657f-40f5-8316-3b31d08ed66e)
