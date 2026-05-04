# Bayesian Multifractal Image Segmentation

Our work introduces an unsupervised Bayesian multifractal segmentation method to model and segment multifractal textures by jointly estimating the multifractal parameters and labels on
images, at the pixel-level.

<p align="center">
  <a href="https://sites.google.com/view/kleonlopez/home">Kareth León</a><sup>1</sup>, 
  <a href="https://abderrahim-halimi.vercel.app/index.html">Abderrahim Halimi</a><sup>2</sup>, 
  <a href="https://perso.tesa.prd.fr/jyt/">Jean-Yves Tourneret</a><sup>1</sup>, and
  <a href="https://www.irit.fr/~Herwig.Wendt/index.html">Herwig Wendt</a><sup>1</sup>

  ><i lang="it"><small><sup>1</sup>IRIT Laboratory, CNRS, INP-Toulouse, UT3, UT2, TéSA, Toulouse, France. <small></i>
  ><i lang="it"><sup>2</sup>Heriot-Watt University, Edinburgh EH14 4AS, UK. </i>  

  <img src="https://github.com/kareth-leon/MFA_Seg/blob/main/image/img_results.png">
  </p>


>In this repository: demo implementation in Matlab for the segmentation of a 2D multifractal random walk (MRW).

## Some Results on Real Data

#### Water detection from SAR images


<blockquote>
<p align="center">
  <img src="https://github.com/kareth-leon/MFA_Seg/blob/main/image/SAR_results.png">
</p>
  >Dataset taken from https://github.com/myeungun/SAR-water-segmentation.
</blockquote>

## How to cite

If you use our code in your research, please cite our work:

```Latex
@article{leon2025bayesian,
  title={Bayesian Multifractal Image Segmentation},
  author={Le{\'o}n-L{\'o}pez, Kareth M and Halimi, Abderrahim and Tourneret, Jean-Yves and Wendt, Herwig},
  journal={IEEE Transactions on Image Processing},
  volume={34},
  pages={8500--8510},
  year={2025},
  publisher={IEEE}
}
```
Paper available in Arxiv also! <a href="https://arxiv.org/pdf/2501.08694">here</a>.



## Acknowledgements
- This work was implemented using the Multifractal Analysis Matlab toolbox <a href="https://www.irit.fr/~Herwig.Wendt/software.html#bayesc2)">"Bayesian univariate and multivariate models and estimators for (c1,c2)"</a> of H. Wendt.
- This work was supported by the Project MUTATION - French National ''ANR JCJC'' Grant - 2019-2023. A. Halimi was supported by the UK Royal Academy of Engineering under the Research Fellowship Scheme (RF/201718/17128). 
