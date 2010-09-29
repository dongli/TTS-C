program TestPolarStereo

    use SphereService
    use MeshManager

    implicit none

    real(8) x0(2), x1(2)

    write(*, "('Input a spherical coordinate:')")
    read(*, *) x0

    x0 = x0/Rad2Deg

    call PolarStereoTransformCoordinate(x0, x1, 1)
    
    write(*, "('Transformed coordinate:')")
    write(*, *) x1

    call PolarStereoTransformCoordinateBack(x1, x0, 1)

    write(*, "('Transformed back coordinate:')")
    write(*, *) x0*Rad2Deg

end program TestPolarStereo

