.PHONY: clean

main: grid.o mpi_interface.o matrix.o init.o netcdf_write.o main.o
	mpifort -o main grid.o mpi_interface.o matrix.o init.o netcdf_write.o main.o \
	-I/usr/include -L/usr/lib/x86_64-gnu-linux/ -lnetcdff -lnetcdf

grid.o: grid.f95
	mpifort -c grid.f95

mpi_interface.o: mpi_interface.f95
	mpifort -c mpi_interface.f95

matrix.o: matrix.f95
	mpifort -c matrix.f95

init.o: init.f95
	mpifort -c init.f95

netcdf_write.o: netcdf_write.f95
	mpifort -c netcdf_write.f95 -I/usr/include -L/usr/lib/x86_64-gnu-linux/ -lnetcdff -lnetcdf

main.o: main.f95
	mpifort -c main.f95

clean:
		rm -f grid.o init.o netcdf_write.o mpi_interface.o matrix.o main.o main grid.mod init.mod netcdf_write.mod mpi_interface.mod matrix.mod

