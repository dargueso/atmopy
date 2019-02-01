
class ObsInfo(object):
    datasets = {
        'CMORPH_CRT':{'path_obs':'/home/dargueso/OBS_DATA/CMORPH/CRT/',
                      'patt_obs':'MC_CMORPH_V1.0_ADJ_8km-',
                      'pr_varname': 'pr',
                      'scale': 3600.,
                      'latvar':'lat',
                      'lonvar':'lon',
                      },

        'CMORPH_RAW':{'path_obs':'/home/dargueso/OBS_DATA/CMORPH/RAW/',
                    'patt_obs':'MC_CMORPH_V1.0_RAW_8km-',
                    'pr_varname': 'pr',
                    'scale': 3600.,
                    'latvar':'lat',
                    'lonvar':'lon',
                    },

        'TRMM':    {'path_obs':'/home/dargueso/OBS_DATA/TRMM/',
                    'patt_obs':'TRMM_',
                    'pr_varname': 'pcp',
                    'scale': 1.,
                    'latvar':'latitude',
                    'lonvar':'longitude',
                    },

        'GPM':    { 'path_obs':'/home/dargueso/OBS_DATA/GPM',
                    'patt_obs':'MC_GPM_0.1_',
                    'pr_varname': 'PR',
                    'scale': 1.,
                    'latvar':'lat',
                    'lonvar':'lon',
                    },
        # 'CHIRPS':{ 'path_obs':'/home/dargueso/OBS_DATA/CHIRPS-2.0',
        #             'patt_obs':'chirps-v2.0.2012-2016_',
        #             'pr_varname': 'precip',
        #             'scale': 3600.,
        #             'latvar':'latitude',
        #             'lonvar':'longitude',
        #             },
                }


    def get_avail_datasets(self):
        """ get the list of available obs datasets
        """

        return self.datasets.keys()

    def get_param(self, obs_name):
        """ get the list of parameters avaiable for a dataset
        """
        return self.datasets[obs_name].keys()

    def get_param_values(self,obs_name,param):
        """ get the value of a given parameters for a dataset
        """

        return self.datasets[obs_name][param]
