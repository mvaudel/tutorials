# Mass spectrometry-based proteomics data

In this tutorial, we navigate mass spectrometry-based proteomics raw data. Please note that there is a vast variety of acquisition modes resulting in differences in the data. Notably, one should distinguish (i) _targeted_ proteomics, where specific compounds are targeted by the mass spectrometer; and (ii) _untargeted_ acquisition, where the mass spectrometer screens in an unbiased fashion. Two acquisition modes can be distinguished in _untargeted_ proteomics: (i) _Data Dependent Acquisition_, DDA, where the mass spectrometer attempts at fragmenting one compound at a time by selecting them based on mass, and (ii) _Data Independent Acquisition_, DIA, where the mass spectrometer fragments all compounds within a given mass range. Here, we will focus on untargeted DDA proteomics data.

Proteomics data consist of (1) MS1 spectra, also called survey scan, the mass spectra of compounds entering the mass spectrometer at a given time point, and (2) MSn spectra, the mass spectra of compounds after n fragmentation steps. In most cases, compounds are fragmented once, resulting in MS1 and MS2 spectra, acquired before and after fragmentation, respectively. As schematized in Figure 1 below, a proteomic data set can be represented in multiple dimensions: the time of elution, m/z and intensities.

![Proteomics Data](images/data.png?raw=true "Proteomics data")
_Figure 1: schematic representation of proteomics data, taken from [(1)](#references)._

Proteomics data consist of multiple mass spectrometry runs, each corresponding to the acquisition of a sample. A sample can be spread across multiple runs, which is referred to as _sample fractionation_. The organization of samples and files depends on the experimental design, and it is important to take the experimental design into account when analysing data.


## Public proteomics data

In a global effort for scientific transparency, scientists worldwide share the proteomics data supporting their conclusion. The ProteomeXchange consortium is a worldwide coordinated effort for sharing proteomics data [(2)](#references).

![ProteomeXchange overview](http://www.proteomexchange.org/px_members.png "ProteomeXchange overview")

:pencil2: Navigate to [ProteomeXchange](http://www.proteomexchange.org), click on `Access Data`.

In this section, you can navigate the data provided by the community and look for data sets of interest.

:pencil2: Open the data set `PXD001819`.

You now see information on the data set, the publication as well as metadata on the project.

:pencil2: Click `PRIDE project URI`.

You are now redirected to PRIDE [(3)](#references), the repository where the data were deposited. You can access more information as well as the original files.

:pencil2: Click `Download Project Files`.

You now see the files as they were provided by the submitter.

:pencil2: Scroll down to `RAW Files` and download `UPS1_12500amol_R1.raw`.


## Format and tools

Mass spectrometers produce raw data in various formats which are not necessarily open. These formats can be converted to mzML, the reference standard for mass spectrometry files [(4)](#references). The reference tool to process raw files is Proteowizzard [(5)](#references), it allows navigating and converting most mass spectrometry formats. Note that the libraries needed to open the files are often provided for Windows only, Mac and Linux users will encounter difficulties working with the raw data.

:pencil2: Install [Proteowizzard](http://proteowizard.sourceforge.net/)


## Raw data

Navigating the raw data will allow you to find specific spectra, and look at the spectra acquired by the instrument. It is the first step of quality control (QC) and can help spotting many mistakes.

:pencil2: Go to the installation folder and open `SeeMS`.

:pencil2: Open the file `UPS1_12500amol_R1.raw` downloaded in the previous section.

You should see the following screen.

![SeeMS Overview](images/seeMS_1.png?raw=true "SeeMS Overview")

The top screen shows the ion chromatogram, it is the sum of MS1 intensities at a given elution time.

![Ion Chromatogram](images/chromatogram.png?raw=true "Ion Chromatogram")

:pencil2: Zoom in a section of the chromatogram

You should see something like this.

![Ion Chromatogram Zoom](images/chromatogram_zoom.png?raw=true "Ion Chromatogram Zoom")

[:thought_balloon:1](../Answers.md#thought_balloon1) _What are the gaps that we observe?_

In the table at the bottom of the screen, you can see all the spectra recorded. You can select spectra and visualize them. In the image below, you can see an MS1 and MS2 specrum selected.

![MS1 MS2](images/MS1_MS2.png?raw=true "MS1 MS2")

[:thought_balloon:2](../Answers.md#thought_balloon2) _In proteomics, what do we expect to find in these spectra?_


## References

(1) [Anatomy and evolution of database search engines-a central component of mass spectrometry based proteomic workflows](https://www.ncbi.nlm.nih.gov/pubmed/28902424)

(2) [ProteomeXchange provides globally coordinated proteomics data submission and dissemination](https://www.ncbi.nlm.nih.gov/pubmed/24727771)

(3) [PRIDE: the proteomics identifications database](https://www.ncbi.nlm.nih.gov/pubmed/16041671)

(4) [mzML-a community standard for mass spectrometry data](https://www.ncbi.nlm.nih.gov/pubmed/20716697)

(5) [A cross-platform toolkit for mass spectrometry and proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23051804)

