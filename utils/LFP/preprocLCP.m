function [lfpDat] = preprocLFP(lfpDat)


cfg                     = [];
cfg.detrend             = 'yes';
cfg.demean              = 'yes';


[lfpDat] = ft_preprocessing( cfg , lfpDat )