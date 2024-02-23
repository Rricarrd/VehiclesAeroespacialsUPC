function fillhdf(filename_template,filename_output,u)
    copyfile(filename_template,filename_output);  
    dset_name = '/NASTRAN/RESULT/NODAL/DISPLACEMENT'; 
    dataset = h5read(filename_output,dset_name);
    dataset.X = u(:,1);
    dataset.Y = u(:,2);
    dataset.Z = u(:,3);
    dataset.RX = u(:,4);
    dataset.RY = u(:,5);
    dataset.RZ = u(:,6);

    fileattrib(filename_output,'+w');
    plist = 'H5P_DEFAULT';
    fid = H5F.open(filename_output,'H5F_ACC_RDWR',plist);
    dset_id = H5D.open(fid,dset_name);
    H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,dataset);
    H5D.close(dset_id);
    H5F.close(fid);
end