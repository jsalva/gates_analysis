
subject_list=['300','301','302','303','304','305','306','307','308','309','310','311','312',
'313','314','315','316','317','318','319','320','321','322','323','325','326','327',
'328','329','330','332','333','334','335','336','337','338','339','340','341','342',#'338',temp exclusion from MCNAB for no behav data
'343','344','345','346','347','348','349','350','351','352','353','354','355','356',
'357','358','359','360','361','362','363','364','365','366','367','368','369','370',
'371','372','373','374','375','376','377','378','379','380','381','382','383','384',
'385','386','387','388','389','390','391','392','393',
'400','401','402','403','404','405','406','407','408','409','410','411','412',
'413','414','415','416',
'500','501','502','503','504','505']#'300','301', excluded because couldnt find data files

info = {}
info['300']= [(['mprag_wrongprotocol1'],'struct'),(['WM2','WM3'],'WM'),(['Nback_199tr1','Nback_199tr2'],'Nback_letters'),([''],'mag'),([''],'phase')]
info['301']= [(['mprag_wrongprotocol2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback_199tr1','Nback_199tr2'],'Nback_letters'),([''],'mag'),([''],'phase')]
#info['400']= [(['mprag1'],'struct'),([na,na,na,na],'WM'),(['NBack1','NBack2'],'Nback_letters'),(['NBack3','NBack4'],'Nback_spatial')]
info['302']=[(['mprag2'],'struct'),(['WM2','WM4'],'WM'),(['Nback_276tr1','Nback1'],'Nback_letters'),([''],'mag'),([''],'phase')]
#only 1 "good" nback run (second run) first had 204 TR; 
info['303']=[(['mprag2'],'struct'),(['WM2','WM3'],'WM'),(['Nback_204tr1','Nback1'],'Nback_letters'),([''],'mag'),([''],'phase')]
info['304']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),([''],'mag'),([''],'phase')]
info['305']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),([''],'mag'),([''],'phase')]
info['306']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['307']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['308']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['309']=[(['mprag2'],'struct'),(['WM1','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#310 Nback run 1 2 outliers for motion
info['310']=[(['mprag2'],'struct'),(['WM1','WM2','WM3'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['311']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['312']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#313 WMrun2 mot outlier; all come close
info['313']=[(['mprag2'],'struct'),(['WM1','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['314']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['315']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['316']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),([''],'mag'),([''],'phase')]
info['317']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['318']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['319']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['320']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1'],'Nback_letters'),([''],'mag'),([''],'phase')]
info['321']=[(['mprag2'],'struct'),(['WM1','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['322']=[(['mprag2'],'struct'),(['WM2','WM3'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['323']=[(['mprag2'],'struct'),(['WM1','WM2','WM3'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['324']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['325']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#326 motion outlier Nback run2
info['326']=[(['mprag2'],'struct'),(['WM5','WM2'],'WM'),(['Nback1'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['327']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['328']=[(['mprag2'],'struct'),(['WM1','WM3'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['329']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['330']=[(['mprag2'],'struct'),(['WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['332']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['333']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['334']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['335']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['336']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['337']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#all nans for 338
info['338']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['339']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['340']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['341']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['342']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['343']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['344']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['345']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['346']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['347']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['348']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#349 motion outlier run2 nback
info['349']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['350']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['351']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#352 run1 and 4 for WM missing onsets (2R and 4ry respectively)
info['352']=[(['mprag2'],'struct'),(['WM2','WM3'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#353 wm Run4 mot outlier; Nback motion outliers (both runs)
info['353']=[(['mprag2'],'struct'),(['WM1','WM2','WM3'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['354']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['355']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['356']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['357']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['358']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['359']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['360']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['361']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['362']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#363 WM run2 mot outlier; Nback motion outliers (both runs)
info['363']=[(['mprag2'],'struct'),(['WM1','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['364']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['365']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['366']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['367']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['368']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['369']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['370']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['371']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['372']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['373']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['374']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['375']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['376']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['377']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['378']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['379']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['380']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['381']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['382']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['383']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['384']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['385']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['386']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['387']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['388']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['389']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['390']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['391']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['392']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['393']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['400']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['401']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['402']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['403']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['404']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['405']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['406']=[(['mprag2'],'struct'),(['WM1','WM2'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['407']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['408']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['409']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['410']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['411']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['412']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['413']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#414 WMrun1 motion outlier
info['414']=[(['mprag2'],'struct'),(['WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['415']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['416']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['500']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['501']=[(['mprag2'],'struct'),(['WM2','WM3'],'WM'),(['Nback1','Nback2'],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['502']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['503']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['504']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['505']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['b500_2']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#501_2 WM run 4 missing 2R onsets
info['b501_2']=[(['mprag2'],'struct'),(['WM2','WM3'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['b502_2']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
info['b503_2']=[(['mprag2'],'struct'),(['WM1','WM2','WM3','WM4'],'WM'),([''],'Nback_letters'),(['fieldmap_func1'],'mag'),(['fieldmap_func2'],'phase')]
#ADD MORE PARTICIPANTS
