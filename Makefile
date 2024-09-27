# Compiler and flags
FC = gfortran
FLAGS = -Wall -O2

# Source and object directories
SRCDIR = modules
OBJDIR = modules

# Object files
OBJ = $(OBJDIR)/m_cells.o $(OBJDIR)/m_field.o $(OBJDIR)/m_fluid.o $(OBJDIR)/streamer.o

# Default target: build the executable
streamer: $(OBJ)
	$(FC) $(FLAGS) -o streamer_1D $(OBJ)

# Compile m_read_in_data module
$(OBJDIR)/m_cells.o: $(SRCDIR)/m_cells.f90
	$(FC) $(FLAGS) -c $(SRCDIR)/m_cells.f90 -o $(OBJDIR)/m_cells.o

# Compile m_field module
$(OBJDIR)/m_field.o: $(SRCDIR)/m_field.f90
	$(FC) $(FLAGS) -c $(SRCDIR)/m_field.f90 -o $(OBJDIR)/m_field.o

# Compile m_fluid module (depends on m_field if it uses it)
$(OBJDIR)/m_fluid.o: $(SRCDIR)/m_fluid.f90
	$(FC) $(FLAGS) -c $(SRCDIR)/m_fluid.f90 -o $(OBJDIR)/m_fluid.o

# Compile the streamer program (depends on m_field and m_fluid)
$(OBJDIR)/streamer.o: $(SRCDIR)/streamer.f90
	$(FC) $(FLAGS) -c $(SRCDIR)/streamer.f90 -o $(OBJDIR)/streamer.o

# Clean target to remove compiled files
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod streamer 
