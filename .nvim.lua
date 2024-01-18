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
        args = { '-l', '${workspaceFolder}/debug/loadfile',
            '-x', '${workspaceFolder}/debug/SingleX_lr.txt',
            '-e', '0.2', '-E', '0', '-E', '0', '-E', '0', '-p' };
        cwd = '${workspaceFolder}/debug';
        terminal = 'integrated';
        console = 'integratedTerminal';
    }
}

