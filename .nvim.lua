local dap = require('dap')
dap.adapters.codelldb = {
    type = 'server';
    port = '${port}';
    executable = {
        command = vim.fn.exepath('codelldb'),
		args = { '--port', '${port}'},
    };
}
dap.configurations.cpp = {
    {
        type = 'codelldb';
        request = "launch";
        name = 'Debug';
        program = '${workspaceFolder}/build/SXCpp';
        args = { '-l', '${workspaceFolder}/debug/Load/Load3e-3_plane.txt',
            '-x', '${workspaceFolder}/debug/SingleX/SingleX_neuler.txt',
            '-s', '0.1', '-o', '0.001', '-b', '0.1',
            '-e', '0.2', '-E', '0', '-E', '0', '-E', '0', '-p' };
        cwd = '${workspaceFolder}/debug/00_vim';
        terminal = 'integrated';
        console = 'integratedTerminal';
    }
}

