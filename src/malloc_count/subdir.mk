################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
./src/malloc_count/malloc_count.c 

OBJS += \
./src/malloc_count/malloc_count.o 

C_DEPS += \
./src/malloc_count/malloc_count.d 


# Each subdirectory must supply rules for building sources it contributes
src/malloc_count/%.o: ./src/malloc_count/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


