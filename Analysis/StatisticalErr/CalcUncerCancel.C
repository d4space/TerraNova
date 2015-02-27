int CalcUncerCancel()
{
  double WpNorm_MomRes[] = {0.352, 0.203, 0.143, 0.705, 0.604, 0.460, 0.355, 0.348, 0.438, 0.510, 0.546};
  double WmNorm_MomRes[] = {0.266, 0.146, 0.088, 0.412, 0.454, 0.495, 0.532, 0.558, 0.579, 0.621, 0.642};

  double WpDiff_MomRes[] = {0.390, 0.210, 0.260, 0.380, 0.660, 0.580, 0.330, 0.350, 0.470, 0.560, 0.460};
  double WmDiff_MomRes[] = {0.150, 0.100, 0.180, 0.240, 0.380, 0.600, 0.240, 0.340, 0.440, 0.520, 0.560};

  double Wp_MomResCancel[11];
  double Wm_MomResCancel[11];

  for(int i(0);i<11;i++)
  {
    Wp_MomResCancel[i] = 1-WpNorm_MomRes[i]/WpDiff_MomRes[i];
    Wm_MomResCancel[i] = 1-WmNorm_MomRes[i]/WmDiff_MomRes[i];
    cout<<"WpMomResBin "<<i+1<<"\t"<<Wp_MomResCancel[i]*100<<"%"<<endl;
    cout<<"WmMomResBin "<<i+1<<"\t"<<Wm_MomResCancel[i]*100<<"%"<<endl;
  }

  double WpNorm_MEtRes[] = {0.055, 0.030, 0.062, 0.098, 0.098, 0.103, 0.122, 0.152, 0.182, 0.205, 0.218};
  double WmNorm_MEtRes[] = {0.033, 0.021, 0.043, 0.061, 0.068, 0.095, 0.146, 0.208, 0.264, 0.307, 0.329};
  double WpDiff_MEtRes[] = {0.060, 0.050, 0.090, 0.140, 0.250, 0.300, 0.120, 0.210, 0.360, 0.550, 0.660};
  double WmDiff_MEtRes[] = {0.040, 0.030, 0.060, 0.090, 0.170, 0.240, 0.150, 0.300, 0.480, 0.650, 0.770};

  double WpNorm_QcdBkg[] = {0.858, 1.364, 1.157, 2.205, 1.264, 1.892, 0.904, 0.910, 0.922, 0.939, 0.957};
  double WmNorm_QcdBkg[] = {0.891, 1.234, 1.309, 1.930, 1.199, 2.292, 1.006, 1.025, 1.040, 1.057, 1.079};
  double WpDiff_QcdBkg[] = {1.190, 1.600, 1.130, 2.260, 1.150, 1.860, 0.690, 0.680, 0.700, 0.730, 0.740};
  double WmDiff_QcdBkg[] = {1.240, 1.410, 1.340, 1.970, 1.080, 2.310, 0.840, 0.840, 0.860, 0.890, 0.910};

  double Wp_QcdBkg[11];
  double Wm_QcdBkg[11];
  for(int i(0);i<11;i++)
  {
    Wp_QcdBkg[i] = 1-WpNorm_QcdBkg[i]/WpDiff_QcdBkg[i];
    Wm_QcdBkg[i] = 1-WmNorm_QcdBkg[i]/WmDiff_QcdBkg[i];
    cout<<"WpQcdBin "<<i+1<<"\t"<<Wp_QcdBkg[i]*100<<"%"<<endl;
    cout<<"WmQcdBin "<<i+1<<"\t"<<Wm_QcdBkg[i]*100<<"%"<<endl;
  }

  double WpNorm_QcdShape[] = {0.260, 0.392, 0.255, 0.410, 0.239, 0.332, 0.381, 0.220, 0.891, 0.947, 0.939};
  double WmNorm_QcdShape[] = {0.127, 0.261, 0.383, 0.224, 0.149, 0.570, 0.272, 0.149, 0.845, 0.895, 0.898};
  double WpDiff_QcdShape[] = {0.370, 0.460, 0.210, 0.390, 0.170, 0.290, 0.350, 0.130, 0.880, 0.940, 0.920};
  double WmDiff_QcdShape[] = {0.130, 0.310, 0.420, 0.210, 0.110, 0.580, 0.260, 0.110, 0.840, 0.900, 0.890};

  double WpNorm_EWK[] = {0.008, 0.015, 0.030, 0.030, 0.063, 0.099, 0.132, 0.159, 0.179, 0.193, 0.199};
  double WmNorm_EWK[] = {0.009, 0.019, 0.041, 0.034, 0.075, 0.108, 0.129, 0.143, 0.153, 0.160, 0.164};
  double WpDiff_EWK[] = {0.170, 0.150, 0.150, 0.190, 0.240, 0.270, 0.340, 0.360, 0.380, 0.390, 0.400};
  double WmDiff_EWK[] = {0.200, 0.190, 0.180, 0.240, 0.300, 0.360, 0.340, 0.350, 0.360, 0.360, 0.370};

  double Wp_EWK[11];
  double Wm_EWK[11];
  for(int i(0);i<11;i++)
  {
    Wp_EWK[i] = 1-WpNorm_EWK[i]/WpDiff_EWK[i];
    Wm_EWK[i] = 1-WmNorm_EWK[i]/WmDiff_EWK[i];
    cout<<"WpEWK "<<i+1<<"\t"<<Wp_EWK[i]*100<<"%"<<endl;
    cout<<"WmEWK "<<i+1<<"\t"<<Wm_EWK[i]*100<<"%"<<endl;
  }


  return 0;
}

