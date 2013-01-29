function returnval = ovrlapr(sCandMaskParadigm, Tthresh)
    spm_defaults                       % I'm not sure if I want this still, though I'm pretty sure it wont hurt.
    returnval = [];
    subjlist = {
    '300','301'
    };
        basedir='/mindhive/gablab/GATES/Analysis/ROI/';  
        resultsdir = 'results/';
        cellTmaps = {'0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008'};
        for i = 1:length(cellTmaps)
            
            filename = [basedir resultsdir 'structLMFG_Ttresh_' num2str(Tthresh) '_' sCandMaskParadigm '_WM' cellTmaps{i} '.txt'];

            [outfid(i).file,fidmsg] = fopen(filename,'w');
            if outfid(i).file == -1
                disp(fidmsg)
            end
            fprintf(outfid(i).file, '%s\t%f\t\n', 'threshold: ', Tthresh);
            fprintf(outfid(i).file, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Subject', 'WM > Thresh','WM < -Thresh', ['Nback_' sCandMaskParadigm ' >Thresh'],['Nback_' sCandMaskParadigm ' < -Thresh'],'ActWM ActNback overlap','ActWM DeactNback overlap', 'DeactWM ActNback overlap', 'DeactWM DeactNback overlap')   
            
        end
    for sub = 1:length(subjlist)
        sSubId = subjlist{sub};
        clear maps
        tdir = ['/mindhive/gablab/GATES/Analysis/WM/l1output/' sSubId '/subj_contrasts/'];
        for i = 1:length(cellTmaps)
            maps(i).str = ['spmT_' cellTmaps{i} '.img'];
            maps(i).tpath = [tdir maps(i).str];
        end
 
        
        sMask = [sSubId '_LMFG_mask.nii'];
       
        maskdir= [basedir sSubId '/'];                 % I think this needs to be changed to the roi dir where these are stored. (look at the saved version I have on my external hard drive)
            % Where the t-map where are looking in is located.
        candTmask = ['/mindhive/gablab/GATES/Analysis/Nback_' sCandMaskParadigm '/l1output/' sSubId '/subj_contrasts/spmT_0001.img'];
        

    %    clTmap = 'spmT_0003.img';          % keeping these numbers the same since these are still the correct t maps 8 = CL_nonerror, 12 = Eng_nonerror
    %    englishTmap = 'spmT_0008.img';     % This is the spmTmap for which the number of supra thresh voxels are being appraised.
        mask = [maskdir sMask]; % This is the name of the mask where we are looking. Was this in the basedir before?
        masktext = sMask;           % For the naming of the file later?
    
%    d = dir([basedir, '3*']);          % Means it will run this for all subjects who start with a 1

%     disp(' ');
%     disp('Available subjects:');
%     for l = 1:length(d)
%         if l <= 10
%             spc = ' ';
%         else
%             spc = [];
%         end
%         switch rem(l,2)
%             case 1, tmp = ['  ', int2str(l), ') ', d(l).name, spc];
%             case 0, disp([tmp, '    ', int2str(l), ') ', d(l).name]);
%         end
%     end
% 
%     done = 0;
%     count = 1;
%     disp('Please select files (-1 for all, 0 to finish):');
%     while ~done & count <= length(d)
%         subject(count) = input('  File number:  ');
%         if subject(count) == -1
%             subject = [1:length(d)];
%             done = 1;
%         elseif subject(count) == 0
%             subject = subject(1:count-1);
%             done = 1;
%         elseif subject(count) < 1 | subject(count) > length(d)
%             disp('    No such selection.');
%         elseif ~isempty(find(subject(count) == subject(1:count-1)))
%             disp('    Selection already made.');
%         else
%             count = count + 1;
%         end
%     end

    %%% get t threshold from user
%    thresh = input('Enter threshold: ');

    %%% name output file


%    for s = 1:length(subject)
%        sub = d(subject(s)).name;

        %%% load up mask
        mask_hdr = spm_vol(mask);
        mask_data = spm_read_vols(mask_hdr);

        % read in contrasts... first cl then english
        for i = 1:length(maps)

            maps(i).hdr = spm_vol(maps(i).tpath);
            maps(i).data = spm_read_vols(maps(i).hdr);
            maps(i).data_inmask = maps(i).data(find(mask_data));
            maps(i).actST = length(find(maps(i).data_inmask>=Tthresh));
            maps(i).deactST = length(find(maps(i).data_inmask<=-Tthresh));
            
            maps(i).maskpath = candTmask;
            maps(i).maskhdr = spm_vol(maps(i).maskpath);
            maps(i).maskdata = spm_read_vols(maps(i).maskhdr);
            maps(i).maskdata_inmask = maps(i).maskdata(find(mask_data));
            maps(i).maskactST = length(find(maps(i).maskdata_inmask>=Tthresh));
            maps(i).maskdeactST = length(find(maps(i).maskdata_inmask<=-Tthresh));
%         cl_str = [tdir clTmap];
%         cl_hdr = spm_vol(cl_str);
%         cl_data = spm_read_vols(cl_hdr);
%         cl_data_inmask = cl_data(find(mask_data));
% 
% 
%         en_str = [tdir englishTmap];
%         en_hdr = spm_vol(en_str);
%         en_data = spm_read_vols(en_hdr);    
%         en_data_inmask = en_data(find(mask_data));


        %%% count super threshold voxels in each tmap and just store in
        %%% variable (clST = created language super threshold)
%         clST = length(find(cl_data_inmask >= thresh));
%         enST = length(find(en_data_inmask >= thresh));

        %%% count number of overlapping superthreshold voxels and store in a
        %%% variable
        maps(i).actactcount = 0;
        maps(i).actdeactcount = 0;
        maps(i).deactactcount = 0;
        maps(i).deactdeactcount=0;
        for j = 1:length(find(mask_data))
            if maps(i).data_inmask(j) >= Tthresh && maps(i).maskdata_inmask(j) >= Tthresh
                maps(i).actactcount = maps(i).actactcount + 1;
            elseif maps(i).data_inmask(j)>=Tthresh && maps(i).maskdata_inmask(j) <=-Tthresh
                maps(i).actdeactcount=maps(i).actdeactcount+1;
            elseif maps(i).data_inmask(j)<=-Tthresh && maps(i).maskdata_inmask(j) >=Tthresh
                maps(i).deactactcount=maps(i).deactactcount + 1;
            elseif maps(i).data_inmask(j)<=-Tthresh && maps(i).maskdata_inmask(j) <=-Tthresh
                maps(i).deactdeactcount=maps(i).deactdeactcount+1;
            end
        end
    
        fprintf(outfid(i).file, '%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', sSubId, maps(i).actST, maps(i).deactST, maps(i).maskactST, maps(i).maskdeactST, maps(i).actactcount, maps(i).actdeactcount, maps(i).deactactcount, maps(i).deactdeactcount)
        disp(sSubId)
        end
    end
    %end
end