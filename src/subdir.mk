################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
./src/CollectionInfo.cpp \
./src/GESAConverter.cpp \
./src/MultiACS.cpp \
./src/Parameters.cpp \
./src/Reader.cpp \
./src/StackedDGenerator.cpp \
./src/Writer.cpp 

OBJS += \
./src/CollectionInfo.o \
./src/GESAConverter.o \
./src/MultiACS.o \
./src/Parameters.o \
./src/Reader.o \
./src/StackedDGenerator.o \
./src/Writer.o 

CPP_DEPS += \
./src/CollectionInfo.d \
./src/GESAConverter.d \
./src/MultiACS.d \
./src/Parameters.d \
./src/Reader.d \
./src/StackedDGenerator.d \
./src/Writer.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ./src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


