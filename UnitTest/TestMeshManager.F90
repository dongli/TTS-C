program TestMeshManager

    use SphereService
    use MeshManager

    implicit none

    integer i
    type(Location) loc

    call MeshManager_Init(128, 60)

    do i = 1, numMeshType
        call meshes(i)%dump
    end do

    call MeshManager_LocationCheck([120.0d0,45.0d0]/Rad2Deg, loc)
    call loc%dump
    
end program TestMeshManager
