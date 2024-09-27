# Compiler and flags
FC = gfortran
FLAGS = -Wall -O2 -g -ffpe-trap=invalid,zero,overflow -pedantic -finit-real=snan

# Source and object directories
SRCDIR = modules
OBJDIR = modules

# Object files
OBJ = $(OBJDIR)/m_cells.o $(OBJDIR)/m_field.o $(OBJDIR)/m_fluid.o $(OBJDIR)/streamer.o $(OBJDIR)/m_types.o

# Default target: build the executable
streamer: $(OBJ)
	$(FC) $(FLAGS) -o streamer_1D $(OBJ)

# Compile m_read_in_data module
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FLAGS) -c $< -o $@

# Clean target to remove compiled files
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod streamer

$(OBJDIR)/streamer.o: $(OBJDIR)/m_cells.o $(OBJDIR)/m_field.o $(OBJDIR)/m_fluid.o $(OBJDIR)/m_types.o
$(OBJDIR)/m_fluid.o: $(OBJDIR)/m_types.o
$(OBJDIR)/m_cells.o: $(OBJDIR)/m_types.o
$(OBJDIR)/m_field.o: $(OBJDIR)/m_types.o



