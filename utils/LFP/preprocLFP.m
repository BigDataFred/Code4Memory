function [lfpDat] = preprocLFP(lfpDat, cfg)

if nargin <2
    cfg                     = [];
    cfg.detrend             = 'yes';
    cfg.demean              = 'yes';
end;

[lfpDat] = ft_preprocessing( cfg , lfpDat );