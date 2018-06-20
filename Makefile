RM := rm -rf

SUBDIRS := \
src/models \
src \
src/hmm \
src/heuristics \
src/core \

CC=g++ 
INCDIR1 = $(CURDIR)/src
INCDIR2 = $(CURDIR)
INC=$(INCDIR1) $(INCDIR2)
INC_PAR=$(foreach d, $(INC), -I$d)
CPPFLAGS=$(INC_PAR) -O3 -c -fmessage-length=0 -std=c++11 -msse2 -mfpmath=sse

-include src/models/subdir.mk
-include src/hmm/subdir.mk
-include src/heuristics/subdir.mk
-include src/core/subdir.mk
-include src/subdir.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif


# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: paHMMgene

# Tool invocations
paHMMgene: $(OBJS) $(USER_OBJS)
	g++  -o "paHMM-gene" $(OBJS) $(USER_OBJS) $(LIBS)

# Other Targets
clean:
	-$(RM) $(C++_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) paHMMgene

.PHONY: all clean dependents
.SECONDARY:
