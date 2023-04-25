function mat_ = load_write_mat(infile,olddname,outfile,newdname)
    mat_     = h5read   (infile,olddname);
   
    h5create  (outfile,newdname,size(mat_));
    h5write   (outfile,newdname,mat_);
    
    if(~strcmp(olddname(end-3:end-2),'ii'))
        neFLR    = h5readatt(infile,olddname,'neFLR');
        h5writeatt(outfile,newdname,'neFLR',neFLR);
    end
    if(~strcmp(olddname(end-3:end-2),'ee'))    
        niFLR    = h5readatt(infile,olddname,'niFLR');
        h5writeatt(outfile,newdname,'niFLR',niFLR);
    end
end
